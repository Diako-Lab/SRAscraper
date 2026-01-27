#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
SRAscraper Rule 1: Metadata Collection

This rule collects metadata from GEO/SRA including:
- Sample information from GEO
- Run accessions from SRA
- CRITICAL: Original filenames from Entrez API (for R1/R2/I1/I2 mapping)

Key output: file_metadata.pkl containing:
- expected_read_types_ordered: ORDERED LIST of read types matching split file order
- original_files: List of original filenames with their read types
"""

import os
import sys
import requests
import time
import json
import pickle
import re
import pandas as pd
from xml.etree import ElementTree as ET
from datetime import datetime
from multiprocessing.pool import ThreadPool
from math import ceil
from time import perf_counter

# Snakemake parameters
metadata_dir = os.path.join(snakemake.params.output_dir, 'metadata')
NCBI_search_txt = snakemake.params.NCBI_search_txt
computing_threads = snakemake.params.computing_threads

pd.options.display.max_colwidth = 10000

os.makedirs(metadata_dir, exist_ok=True)
os.chdir(metadata_dir)


# ============================================================================
# NCBI SRA ENTREZ API - PRIMARY SOURCE FOR ORIGINAL FILENAMES
# ============================================================================

def get_sra_metadata_entrez(run_accession, max_retries=3):
    """
    Query NCBI SRA Entrez API and extract:
    - Original submitted filenames (SRAFile with supertype="Original")
    - Library information (strategy, layout, source)
    - Platform information
    
    CRITICAL: Returns expected_read_types_ordered as an ORDERED LIST
    that matches the order of split files from fastq-dump.
    
    The order of original files in Entrez XML corresponds to READ_INDEX order,
    which determines the _1, _2, _3, _4 split file numbering.
    """
    url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
    params = {
        'db': 'sra',
        'id': run_accession,
        'rettype': 'full',
        'retmode': 'xml'
    }
    
    for attempt in range(max_retries):
        try:
            response = requests.get(url, params=params, timeout=30)
            
            if response.status_code != 200:
                time.sleep(1)
                continue
            
            root = ET.fromstring(response.text)
            
            result = {
                'run_accession': run_accession,
                'library_strategy': None,
                'library_layout': None,
                'library_source': None,
                'platform': None,
                'instrument_model': None,
                'original_files': [],           # List of file info dicts (IN ORDER)
                'expected_read_types_ordered': [],  # ORDERED list: [I1, I2, R1, R2]
                'file_types_present': set(),    # Set for quick lookups
                'is_10x': False,
                'is_bam': False,
                'num_expected_files': 0
            }
            
            # Parse library info
            for lib in root.iter('LIBRARY_DESCRIPTOR'):
                strategy = lib.find('LIBRARY_STRATEGY')
                source = lib.find('LIBRARY_SOURCE')
                layout = lib.find('LIBRARY_LAYOUT')
                
                if strategy is not None:
                    result['library_strategy'] = strategy.text
                if source is not None:
                    result['library_source'] = source.text
                if layout is not None:
                    if layout.find('PAIRED') is not None:
                        result['library_layout'] = 'PAIRED'
                    elif layout.find('SINGLE') is not None:
                        result['library_layout'] = 'SINGLE'
            
            # Parse platform
            for platform in root.iter('PLATFORM'):
                for child in platform:
                    result['platform'] = child.tag
                    instrument = child.find('INSTRUMENT_MODEL')
                    if instrument is not None:
                        result['instrument_model'] = instrument.text
                    break
            
            # Parse SRAFile elements - CRITICAL FOR ORIGINAL FILENAMES
            # The ORDER of files here matches READ_INDEX order (split file order)
            for sra_file in root.iter('SRAFile'):
                supertype = sra_file.get('supertype')
                filename = sra_file.get('filename')
                
                # Only care about Original files (submitted by user)
                if supertype == 'Original' and filename:
                    # Check for S3 URLs in Alternatives child element
                    s3_url = None
                    for alt in sra_file.iter('Alternatives'):
                        alt_url = alt.get('url', '')
                        if alt_url.startswith('s3://'):
                            s3_url = alt_url
                            break
                    
                    file_info = {
                        'filename': filename,
                        'semantic_name': sra_file.get('semantic_name'),
                        'size': sra_file.get('size'),
                        'md5': sra_file.get('md5'),
                        'url': sra_file.get('url'),
                        's3_url': s3_url,
                        'cluster': sra_file.get('cluster'),
                        'read_type': extract_read_type(filename),
                        'is_bam': filename.lower().endswith('.bam'),
                        'is_fastq': is_fastq_filename(filename)
                    }
                    
                    result['original_files'].append(file_info)
                    
                    # Build ORDERED list of read types
                    if file_info['is_fastq'] and file_info['read_type']:
                        result['expected_read_types_ordered'].append(file_info['read_type'])
                        result['file_types_present'].add(file_info['read_type'])
                    
                    if file_info['is_bam']:
                        result['is_bam'] = True
            
            # Number of expected FASTQ files
            result['num_expected_files'] = len([f for f in result['original_files'] if f['is_fastq']])
            
            # Convert set to list for JSON serialization
            result['file_types_present'] = list(result['file_types_present'])
            
            # Detect 10x single-cell
            result['is_10x'] = detect_10x(result)
            
            return result
            
        except Exception as e:
            print(f"    Attempt {attempt + 1} failed: {e}")
            time.sleep(1)
    
    return {
        'run_accession': run_accession,
        'error': 'Failed to fetch metadata',
        'expected_read_types_ordered': [],
        'original_files': []
    }


def extract_read_type(filename):
    """
    Extract R1/R2/I1/I2 from filename.
    Returns None if not found.
    """
    if not filename:
        return None
    
    name_upper = filename.upper()
    
    # Order matters - check more specific patterns first
    patterns = [
        (r'[_\.\-]R1[_\.\-]', 'R1'),
        (r'[_\.\-]R2[_\.\-]', 'R2'),
        (r'[_\.\-]I1[_\.\-]', 'I1'),
        (r'[_\.\-]I2[_\.\-]', 'I2'),
        (r'_R1\.', 'R1'),
        (r'_R2\.', 'R2'),
        (r'_I1\.', 'I1'),
        (r'_I2\.', 'I2'),
        (r'_R1_', 'R1'),
        (r'_R2_', 'R2'),
        (r'_I1_', 'I1'),
        (r'_I2_', 'I2'),
    ]
    
    for pattern, read_type in patterns:
        if re.search(pattern, name_upper):
            return read_type
    
    return None


def is_fastq_filename(filename):
    """Check if filename indicates a FASTQ file"""
    if not filename:
        return False
    name_lower = filename.lower()
    return any(name_lower.endswith(ext) for ext in ['.fastq', '.fq', '.fastq.gz', '.fq.gz'])


def detect_10x(result):
    """
    Detect if this is likely 10x single-cell data.
    """
    # Check library source
    source = (result.get('library_source') or '').lower()
    if 'single cell' in source or 'transcriptomic single cell' in source:
        return True
    
    # Check instrument model
    instrument = (result.get('instrument_model') or '').lower()
    if 'chromium' in instrument or '10x' in instrument:
        return True
    
    # Check file pattern - 10x typically has R1 + R2 + optionally I1/I2
    file_types = set(result.get('file_types_present', []))
    if 'R1' in file_types and 'R2' in file_types:
        if 'I1' in file_types or 'I2' in file_types:
            return True
    
    return False


# ============================================================================
# ORIGINAL GEO/SRA METADATA LOGIC
# ============================================================================

X = pd.read_fwf(NCBI_search_txt, header=None)
ftp_links = X[X[0].str.startswith('FTP download')]
ftp_links = ftp_links.reset_index(drop=True)
ftp_links = ftp_links.loc[0:, 0]

ftp_list = []
for index, row in ftp_links.items():
    try:
        tmp = row.split(') ', -1)[1]
        ftp_list.append(tmp)
    except:
        pass


def url_ok(foo_url):
    """Check if GEO FTP URL is accessible"""
    transformed_url = 'https' + foo_url[3:] + 'soft/'
    max_retries = 10
    retry_delay = 1
    
    for attempt in range(max_retries):
        try:
            response = requests.head(transformed_url, timeout=10)
            if response.status_code == 200:
                return True
        except Exception as e:
            pass
        if attempt < max_retries - 1:
            time.sleep(retry_delay)
    return False


# Multiprocessed URL checking
ftp_list_input = ftp_list
N = len(ftp_list_input)

with ThreadPool(computing_threads) as pool:
    chunksize = ceil(len(ftp_list_input) / computing_threads)
    start = perf_counter()
    results = list(pool.map_async(url_ok, ftp_list_input, chunksize=chunksize).get())
    end = perf_counter()
    print(f'URL checking duration: {end-start:.4f}s')


def countX(lst, x):
    return sum(1 for ele in lst if ele == x)


print(f'''Of the initial {len(ftp_links)} datasets from NCBI search,
      {len(ftp_list)} point to GEO FTP pages.
      {countX(results, True)} returned 200 GET statuses.''')

from itertools import compress
ftp_list = list(compress(ftp_list, results))


def download_with_retry(url, max_retries=10, retry_delay=1):
    import wget
    last_exception = None
    for attempt in range(max_retries):
        try:
            filename = wget.download(url)
            return filename
        except Exception as e:
            last_exception = e
            if attempt < max_retries - 1:
                time.sleep(retry_delay)
    raise last_exception if last_exception else Exception("Download failed")


# ============================================================================
# MAIN METADATA COLLECTION
# ============================================================================

import wget
import GEOparse
import subprocess

subprocess.run(["pysradb", "--version"])

from pysradb.search import GeoSearch
from pysradb.sraweb import SRAweb
sradb = SRAweb()

gse_dict = {}
file_metadata = {}

for url in ftp_list:
    print(f"\n{'='*60}")
    print(f"Processing: {url}")
    print(f"{'='*60}")
    
    try:
        wget_this = url + 'soft/' + url[:-1].split('/', -1)[-1] + '_family.soft.gz'
        filename = download_with_retry(wget_this)
        gse = GEOparse.get_GEO(filepath=filename)
        gse_df = gse.phenotype_data
        gse_id = url[:-1].split('/', -1)[-1]
        
        try:
            columns_to_select = [
                'title', 'geo_accession', 'extract_protocol_ch1',
                'description', 'data_processing', 'platform_id',
                'contact_name', 'contact_department', 'contact_institute',
                'contact_address', 'contact_city', 'contact_state',
                'contact_zip/postal_code', 'contact_country'
            ]
            gse_df_filtered = gse_df.drop(
                columns=gse_df.columns.difference(columns_to_select),
                errors='ignore'
            )
            
            instance = GeoSearch(2, 200, geo_query=gse_id)
            instance.search()
            instance_df = instance.get_df()
            
            srp_dfs = []
            for srp_id in instance_df.study_accession.unique():
                srp_df = sradb.sra_metadata(srp_id, detailed=True)
                srp_dfs.append(srp_df)
            
            if len(srp_dfs) > 1:
                final_srp_df = pd.concat(srp_dfs, ignore_index=True)
            else:
                final_srp_df = srp_dfs[0]
            
            instance_df.drop(
                columns=instance_df.columns.difference(['run_1_accession', 'sample_alias']),
                inplace=True
            )
            instance_df = instance_df.rename(columns={
                'run_1_accession': 'run_accession',
                'sample_alias': 'geo_accession'
            })
            
            final_srp_df = pd.merge(final_srp_df, instance_df, on='run_accession', how='outer')
            
            merged_df = gse_df_filtered.merge(
                final_srp_df,
                right_on='geo_accession',
                left_index=True,
                how='outer'
            )
            merged_df.drop(columns=['geo_accession_x', 'geo_accession_y'], inplace=True, errors='ignore')
            merged_df.dropna(subset=['run_accession'], inplace=True)
            
            gse_dict[gse_id] = merged_df
            
            # ================================================================
            # Get original filenames from SRA Entrez XML
            # ================================================================
            print(f"\nFetching original filenames for {len(merged_df)} runs...")
            
            for run_acc in merged_df['run_accession']:
                print(f"  {run_acc}...", end=" ")
                meta = get_sra_metadata_entrez(run_acc)
                file_metadata[run_acc] = meta
                
                # Summary output
                if meta.get('error'):
                    print(f"ERROR")
                elif meta.get('original_files'):
                    ordered_types = meta.get('expected_read_types_ordered', [])
                    print(f"Found {len(meta['original_files'])} files, ordered types: {ordered_types}")
                else:
                    print(f"No original files")
                
                time.sleep(0.3)  # Rate limiting
            
        except Exception as e:
            print(f"Error processing {gse_id}: {e}")
            gse_df = gse.phenotype_data
            gse_dict[gse_id] = gse_df
            
    except Exception as e:
        print(f"Failed to process {url}: {e}")


# ============================================================================
# SUMMARY AND SAVE
# ============================================================================

num_rows = [len(gse_dict[key]) for key in gse_dict.keys()]

print(f'''\n{'='*60}
METADATA SUMMARY
{'='*60}
Total samples: {sum(num_rows)}
Runs with metadata: {len(file_metadata)}
''')

# Categorize by what we found
has_ordered_types = 0
no_patterns = 0
is_bam = 0
is_10x = 0

for run, meta in file_metadata.items():
    if meta.get('is_bam'):
        is_bam += 1
    if meta.get('is_10x'):
        is_10x += 1
    if meta.get('expected_read_types_ordered'):
        has_ordered_types += 1
    else:
        no_patterns += 1

print(f"Has ordered read types (can map by position): {has_ordered_types}")
print(f"No R1/R2 patterns (needs content analysis): {no_patterns}")
print(f"BAM files: {is_bam}")
print(f"Likely 10x single-cell: {is_10x}")

# Show example of ordered types
print(f"\nExample expected_read_types_ordered:")
for run, meta in list(file_metadata.items())[:3]:
    ordered = meta.get('expected_read_types_ordered', [])
    print(f"  {run}: {ordered}")

# Save
with open('dictionary_file.pkl', 'wb') as f:
    pickle.dump(gse_dict, f)
    print('\nSaved: dictionary_file.pkl')

with open('file_metadata.pkl', 'wb') as f:
    pickle.dump(file_metadata, f)
    print('Saved: file_metadata.pkl')

with open('file_metadata.json', 'w') as f:
    json.dump(file_metadata, f, indent=2, default=str)
    print('Saved: file_metadata.json')

print(f"\n{'='*60}")
print("Metadata collection complete!")
print(f"{'='*60}")

sys.exit(0)
