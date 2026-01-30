#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
SRAscraper Rule 2: Download and Rename

This rule downloads SRA data and renames files to Cell Ranger format.

STRATEGY:
1. Load expected_read_types_ordered from Rule 1 (Entrez API)
2. Check if files already exist (skip if previously downloaded)
3. Download with prefetch + parallel-fastq-dump
4. Map split files BY POSITION: _1→first type, _2→second type, etc.
5. Use content analysis for VALIDATION (warn on mismatch)
6. Rename to Cell Ranger format
7. Update metadata dictionary with only successful samples

The key insight: Entrez original file order matches split file order!
- _1.fastq.gz = first file in Entrez list
- _2.fastq.gz = second file in Entrez list
- etc.

Cell Ranger format: {accession}_S1_L001_{R1|R2|I1|I2}_001.fastq.gz
"""

import os
import sys
import gzip
import subprocess
import pickle
import json
import re
import shutil
from collections import Counter
from datetime import datetime
from pathlib import Path
import copy

# Snakemake parameters
output_dir = snakemake.params.output_dir
metadata_dir = os.path.join(snakemake.params.output_dir, 'metadata')
computing_threads = snakemake.params.computing_threads

os.chdir(metadata_dir)


# ============================================================================
# CONFIGURATION
# ============================================================================

# Read length thresholds for 10x content classification (for validation)
BARCODE_UMI_LENGTHS = {24, 26, 28}  # R1 in 10x
SAMPLE_INDEX_MAX_LENGTH = 12        # I1/I2 in 10x
CDNA_MIN_LENGTH = 50                # R2 in 10x
SAMPLE_SIZE = 500


# ============================================================================
# SKIP DETECTION - CHECK IF ALREADY DOWNLOADED
# ============================================================================

def check_already_downloaded(sample_dir, accession, expected_read_types=None):
    """
    Check if Cell Ranger formatted files already exist for this sample.
    
    Returns:
        dict with:
            - 'complete': True if all expected files exist
            - 'partial': True if some files exist
            - 'existing_files': list of existing Cell Ranger formatted files
            - 'missing_types': list of missing read types
    """
    result = {
        'complete': False,
        'partial': False,
        'existing_files': [],
        'missing_types': [],
        'found_types': []
    }
    
    if not os.path.exists(sample_dir):
        return result
    
    # Look for Cell Ranger formatted files
    cellranger_pattern = re.compile(rf'{accession}_S1_L001_(R1|R2|I1|I2)_001\.fastq\.gz')
    
    for f in os.listdir(sample_dir):
        match = cellranger_pattern.match(f)
        if match:
            read_type = match.group(1)
            filepath = os.path.join(sample_dir, f)
            # Verify file is not empty
            if os.path.getsize(filepath) > 0:
                result['existing_files'].append(f)
                result['found_types'].append(read_type)
    
    # Determine completeness
    if expected_read_types:
        expected_set = set(expected_read_types)
        found_set = set(result['found_types'])
        result['missing_types'] = list(expected_set - found_set)
        result['complete'] = (expected_set == found_set) and len(found_set) > 0
        result['partial'] = len(found_set) > 0 and not result['complete']
    else:
        # No expected types known - consider complete if we have R1 and R2 at minimum
        has_r1 = 'R1' in result['found_types']
        has_r2 = 'R2' in result['found_types']
        result['complete'] = has_r1 and has_r2
        result['partial'] = len(result['found_types']) > 0 and not result['complete']
    
    return result


# ============================================================================
# CONTENT ANALYSIS FUNCTIONS (for validation)
# ============================================================================

def sample_fastq_lengths(filepath, n_reads=SAMPLE_SIZE):
    """Sample reads from FASTQ and get length distribution"""
    sequences = []
    
    try:
        opener = gzip.open if filepath.endswith('.gz') else open
        with opener(filepath, 'rt') as f:
            line_count = 0
            for line in f:
                line_count += 1
                if line_count % 4 == 2:
                    sequences.append(line.strip())
                if len(sequences) >= n_reads:
                    break
    except Exception as e:
        return {'error': str(e), 'dominant_length': None}
    
    if not sequences:
        return {'error': 'No sequences', 'dominant_length': None}
    
    lengths = Counter(len(s) for s in sequences)
    dominant = max(lengths.keys(), key=lambda x: lengths[x])
    
    return {
        'dominant_length': dominant,
        'length_distribution': dict(lengths),
        'n_sequences': len(sequences)
    }


def classify_read_by_length(length):
    """Classify 10x read type by length (for validation)"""
    if length is None:
        return 'UNKNOWN'
    if length in BARCODE_UMI_LENGTHS:
        return 'BARCODE_UMI'  # Typically R1 in 10x (28bp)
    if length <= SAMPLE_INDEX_MAX_LENGTH:
        return 'SAMPLE_INDEX'  # I1/I2 (8-10bp)
    if length >= CDNA_MIN_LENGTH:
        return 'CDNA'  # Typically R2 in 10x (90-150bp)
    return 'UNKNOWN'


def get_expected_content_type(read_type):
    """Map read type to expected content type"""
    if read_type == 'R1':
        return ['BARCODE_UMI', 'CDNA']  # R1 can be either (28bp or 150bp)
    elif read_type == 'R2':
        return ['CDNA']  # R2 is always cDNA (150bp)
    elif read_type in ['I1', 'I2']:
        return ['SAMPLE_INDEX']  # Index reads (8-10bp)
    return ['UNKNOWN']


def analyze_split_files(sample_dir, accession):
    """
    Analyze all split files from fastq-dump.
    Returns list sorted by split number.
    """
    results = []
    
    for f in os.listdir(sample_dir):
        if f.startswith(accession) and f.endswith('.fastq.gz'):
            match = re.match(rf'{accession}_(\d+)\.fastq\.gz', f)
            if match:
                split_num = int(match.group(1))
                filepath = os.path.join(sample_dir, f)
                
                length_info = sample_fastq_lengths(filepath)
                content_type = classify_read_by_length(length_info.get('dominant_length'))
                
                results.append({
                    'filename': f,
                    'split_number': split_num,
                    'filepath': filepath,
                    'dominant_length': length_info.get('dominant_length'),
                    'content_type': content_type
                })
    
    # CRITICAL: Sort by split number to match Entrez order
    results.sort(key=lambda x: x['split_number'])
    
    return results


# ============================================================================
# POSITIONAL MAPPING (PRIMARY METHOD)
# ============================================================================

def map_by_position(split_files, expected_read_types_ordered):
    """
    Map split files to read types BY POSITION.
    
    The Entrez original file order matches the split file order:
    - _1.fastq.gz = first type in expected_read_types_ordered
    - _2.fastq.gz = second type
    - etc.
    
    Returns: dict {filename: read_type}
    """
    mapping = {}
    
    for i, file_info in enumerate(split_files):
        if i < len(expected_read_types_ordered):
            mapping[file_info['filename']] = expected_read_types_ordered[i]
        else:
            # More files than expected types - use generic naming
            mapping[file_info['filename']] = f'R{i+1}'
    
    return mapping


def validate_mapping(mapping, split_files):
    """
    Validate that content matches expected read types.
    Returns list of warnings (empty if all OK).
    """
    warnings = []
    
    # Create lookup for file info
    file_info_lookup = {f['filename']: f for f in split_files}
    
    for filename, read_type in mapping.items():
        info = file_info_lookup.get(filename, {})
        content_type = info.get('content_type')
        expected_content = get_expected_content_type(read_type)
        
        if content_type and content_type not in expected_content:
            length = info.get('dominant_length')
            warnings.append(
                f"{filename} mapped to {read_type}, but content is {content_type} ({length}bp). "
                f"Expected: {expected_content}"
            )
    
    return warnings


# ============================================================================
# FALLBACK: CONTENT-BASED MAPPING (when no metadata)
# ============================================================================

def map_by_content_analysis(split_files):
    """
    Map split files by content analysis when no metadata available.
    Uses read length to classify files.
    """
    mapping = {}
    
    index_count = 0
    has_r1 = False
    has_r2 = False
    
    for f in split_files:
        content_type = f['content_type']
        filename = f['filename']
        
        if content_type == 'BARCODE_UMI' and not has_r1:
            mapping[filename] = 'R1'
            has_r1 = True
        elif content_type == 'CDNA':
            if not has_r2:
                mapping[filename] = 'R2'
                has_r2 = True
            elif not has_r1:
                # Both R1 and R2 are cDNA (no barcode file)
                mapping[filename] = 'R1'
                has_r1 = True
            else:
                mapping[filename] = f'R{f["split_number"]}'
        elif content_type == 'SAMPLE_INDEX':
            index_count += 1
            mapping[filename] = f'I{index_count}'
        else:
            mapping[filename] = f'R{f["split_number"]}'
    
    return mapping


# ============================================================================
# RENAME FILES TO CELL RANGER FORMAT
# ============================================================================

def rename_files_to_cellranger(sample_dir, accession, mapping, split_files):
    """
    Rename files to Cell Ranger format.
    Cell Ranger format: {accession}_S1_L001_{R1|R2|I1|I2}_001.fastq.gz
    """
    renames = []
    file_info_lookup = {f['filename']: f for f in split_files}
    
    for old_name, read_type in mapping.items():
        new_name = f"{accession}_S1_L001_{read_type}_001.fastq.gz"
        old_path = os.path.join(sample_dir, old_name)
        new_path = os.path.join(sample_dir, new_name)
        
        info = file_info_lookup.get(old_name, {})
        
        rename_record = {
            'original': old_name,
            'new': new_name,
            'read_type': read_type,
            'content_type': info.get('content_type'),
            'dominant_length': info.get('dominant_length'),
            'old_path': old_path,
            'new_path': new_path,
            'status': None
        }
        
        try:
            if os.path.exists(old_path):
                if os.path.exists(new_path) and old_path != new_path:
                    os.remove(new_path)
                os.rename(old_path, new_path)
                rename_record['status'] = 'SUCCESS'
            else:
                rename_record['status'] = 'FILE_NOT_FOUND'
        except Exception as e:
            rename_record['status'] = 'ERROR'
            rename_record['error'] = str(e)
        
        renames.append(rename_record)
    
    return renames


# ============================================================================
# DOWNLOAD FUNCTIONS
# ============================================================================

def download_with_prefetch(accession, sample_dir, threads):
    """Download using prefetch + parallel-fastq-dump"""
    result = {'success': False, 'files': [], 'error': None}
    
    parent_dir = os.path.dirname(sample_dir)
    
    # PREFETCH
    print(f"    Prefetching {accession}...")
    prefetch_result = subprocess.run(
        ["prefetch", accession, "-O", parent_dir, "--max-size", "u"],
        capture_output=True, text=True, timeout=3600
    )
    
    if prefetch_result.returncode != 0:
        result['error'] = f"Prefetch failed: {prefetch_result.stderr[:200]}"
        return result
    
    # Find SRA file
    sra_file = os.path.join(sample_dir, f"{accession}.sra")
    if not os.path.exists(sra_file):
        sra_file = os.path.join(sample_dir, accession, f"{accession}.sra")
    if not os.path.exists(sra_file):
        sra_file = os.path.join(parent_dir, accession, f"{accession}.sra")
    
    if not os.path.exists(sra_file):
        result['error'] = "SRA file not found after prefetch"
        return result
    
    print(f"    Found: {sra_file}")
    
    # PARALLEL-FASTQ-DUMP
    print(f"    Extracting with parallel-fastq-dump (threads={threads})...")
    extract_result = subprocess.run(
        ["parallel-fastq-dump", "-s", sra_file,
         "--threads", str(threads),
         "--tmpdir", sample_dir,
         "--outdir", sample_dir,
         "--split-spot", "--split-files", "--gzip"],
        capture_output=True, text=True, timeout=14400
    )
    
    if extract_result.returncode != 0:
        print(f"    Warning: {extract_result.stderr[:200]}")
    
    # List extracted files
    for f in sorted(os.listdir(sample_dir)):
        if f.endswith('.fastq.gz') and f.startswith(accession):
            result['files'].append(f)
    
    # Cleanup SRA file
    try:
        if os.path.exists(sra_file):
            os.remove(sra_file)
        nested_dir = os.path.join(sample_dir, accession)
        if os.path.exists(nested_dir):
            shutil.rmtree(nested_dir)
    except:
        pass
    
    result['success'] = len(result['files']) > 0
    print(f"    Downloaded: {result['files']}")
    
    return result


# ============================================================================
# MAIN PROCESSING
# ============================================================================

# Load metadata from Rule 1
with open('dictionary_file.pkl', 'rb') as f:
    gse_dict_raw = pickle.load(f)
    print('Loaded: dictionary_file.pkl')

# Validate gse_dict - ensure all entries have 'run_accession' column
gse_dict = {}
for key, df in gse_dict_raw.items():
    if 'run_accession' not in df.columns:
        print(f"WARNING: Skipping {key} - no 'run_accession' column found")
        print(f"  Available columns: {list(df.columns)}")
        continue
    if len(df) == 0:
        print(f"WARNING: Skipping {key} - empty dataframe")
        continue
    # Drop any rows with missing run_accession
    df_clean = df.dropna(subset=['run_accession']).copy()
    if len(df_clean) == 0:
        print(f"WARNING: Skipping {key} - no valid run_accession values")
        continue
    gse_dict[key] = df_clean
    print(f"  {key}: {len(df_clean)} samples")

if len(gse_dict) == 0:
    print("ERROR: No valid datasets found with run_accession column!")
    print("Please check your metadata collection step.")
    # Create QC directory and placeholder report to allow pipeline to complete
    qc_dir = os.path.join(output_dir, 'QC')
    os.makedirs(qc_dir, exist_ok=True)
    multiqc_report = os.path.join(qc_dir, 'multiqc_report.html')
    with open(multiqc_report, 'w') as f:
        f.write(f"""<!DOCTYPE html>
<html>
<head><title>MultiQC Report - No Data</title></head>
<body>
<h1>No Valid Data Found</h1>
<p>Timestamp: {datetime.now().isoformat()}</p>
<p>The metadata collection step did not find any datasets with valid run_accession values.</p>
<p>Please check your NCBI search text file and ensure the GEO datasets have linked SRA data.</p>
</body>
</html>""")
    print(f"Created placeholder report: {multiqc_report}")
    sys.exit(0)

print(f"Total valid datasets: {len(gse_dict)}")

file_metadata = {}
if os.path.exists('file_metadata.pkl'):
    with open('file_metadata.pkl', 'rb') as f:
        file_metadata = pickle.load(f)
        print(f'Loaded: file_metadata.pkl ({len(file_metadata)} runs)')
else:
    print('WARNING: No file_metadata.pkl - will use content analysis only')

# Create directories
log_dir = os.path.join(output_dir, 'logs')
qc_dir = os.path.join(output_dir, 'QC')
os.makedirs(log_dir, exist_ok=True)
os.makedirs(qc_dir, exist_ok=True)

# Processing log
processing_log = {
    'timestamp': datetime.now().isoformat(),
    'total_samples': 0,
    'samples': [],
    'summary': {
        'success': 0,
        'partial': 0,
        'failed': 0,
        'skipped_existing': 0,
        'used_positional_mapping': 0,
        'used_content_analysis': 0,
        'validation_warnings': 0
    }
}

# Track successful accessions for each BioProject
successful_accessions = {}  # {bioproject: [accession, ...]}


for key in gse_dict.keys():
    print(f"\n{'='*70}")
    print(f"Processing BioProject: {key}")
    print(f"{'='*70}")
    
    successful_accessions[key] = []
    
    for accession in gse_dict[key]['run_accession']:
        processing_log['total_samples'] += 1
        
        sample_log = {
            'accession': accession,
            'bioproject': key,
            'status': None,
            'method': None,
            'expected_types': [],
            'mapping': {},
            'renames': [],
            'validation_warnings': [],
            'issues': []
        }
        
        print(f"\n--- {accession} ---")
        
        sample_dir = os.path.join(output_dir, 'fastq', key, accession)
        os.makedirs(sample_dir, exist_ok=True)
        
        # Get metadata for this run
        meta = file_metadata.get(accession, {})
        expected_types = meta.get('expected_read_types_ordered', [])
        is_10x = meta.get('is_10x', False)
        
        sample_log['expected_types'] = expected_types
        
        print(f"  Expected types (ordered): {expected_types}")
        print(f"  Is 10x: {is_10x}")
        
        # ================================================================
        # STEP 0: Check if already downloaded
        # ================================================================
        download_check = check_already_downloaded(sample_dir, accession, expected_types)
        
        if download_check['complete']:
            print(f"  ✓ SKIPPING - Already downloaded: {download_check['existing_files']}")
            sample_log['status'] = 'SKIPPED_EXISTING'
            sample_log['method'] = 'PREVIOUSLY_DOWNLOADED'
            sample_log['issues'].append(f"Using existing files: {download_check['existing_files']}")
            processing_log['samples'].append(sample_log)
            processing_log['summary']['skipped_existing'] += 1
            processing_log['summary']['success'] += 1
            successful_accessions[key].append(accession)
            continue
        
        if download_check['partial']:
            print(f"  ⚠ Partial download found: {download_check['existing_files']}")
            print(f"    Missing types: {download_check['missing_types']}")
            print(f"    Re-downloading...")
            # Remove partial files and start fresh
            for f in download_check['existing_files']:
                try:
                    os.remove(os.path.join(sample_dir, f))
                except:
                    pass
        
        # ================================================================
        # STEP 1: Download
        # ================================================================
        download_result = download_with_prefetch(accession, sample_dir, computing_threads)
        
        if not download_result['success']:
            print(f"  ERROR: {download_result.get('error')}")
            sample_log['status'] = 'DOWNLOAD_FAILED'
            sample_log['issues'].append(download_result.get('error'))
            processing_log['samples'].append(sample_log)
            processing_log['summary']['failed'] += 1
            continue
        
        # ================================================================
        # STEP 2: Analyze content (for validation)
        # ================================================================
        print(f"  Analyzing file content...")
        split_files = analyze_split_files(sample_dir, accession)
        
        for f in split_files:
            print(f"    {f['filename']}: {f['dominant_length']}bp → {f['content_type']}")
        
        # ================================================================
        # STEP 3: Create mapping
        # ================================================================
        if expected_types and len(expected_types) == len(split_files):
            # POSITIONAL MAPPING (primary method)
            print(f"  Using POSITIONAL mapping (Entrez order)...")
            mapping = map_by_position(split_files, expected_types)
            sample_log['method'] = 'POSITIONAL'
            processing_log['summary']['used_positional_mapping'] += 1
            
            # Validate
            warnings = validate_mapping(mapping, split_files)
            if warnings:
                sample_log['validation_warnings'] = warnings
                processing_log['summary']['validation_warnings'] += 1
                for w in warnings:
                    print(f"    ⚠ Warning: {w}")
        
        elif expected_types and len(expected_types) != len(split_files):
            # File count mismatch - use content analysis
            print(f"  File count mismatch: expected {len(expected_types)}, got {len(split_files)}")
            print(f"  Using CONTENT-BASED mapping...")
            mapping = map_by_content_analysis(split_files)
            sample_log['method'] = 'CONTENT_ANALYSIS'
            sample_log['issues'].append(f"File count mismatch: expected {len(expected_types)}, got {len(split_files)}")
            processing_log['summary']['used_content_analysis'] += 1
        
        else:
            # No metadata - use content analysis
            print(f"  No metadata - using CONTENT-BASED mapping...")
            mapping = map_by_content_analysis(split_files)
            sample_log['method'] = 'CONTENT_ANALYSIS'
            sample_log['issues'].append("No expected_read_types_ordered in metadata")
            processing_log['summary']['used_content_analysis'] += 1
        
        sample_log['mapping'] = mapping
        
        print(f"  Mapping:")
        for i, (old_name, read_type) in enumerate(mapping.items()):
            info = next((f for f in split_files if f['filename'] == old_name), {})
            expected = expected_types[i] if i < len(expected_types) else '-'
            print(f"    _{i+1} ({info.get('dominant_length')}bp) → {read_type} (expected: {expected})")
        
        # ================================================================
        # STEP 4: Rename files
        # ================================================================
        print(f"  Renaming files...")
        renames = rename_files_to_cellranger(sample_dir, accession, mapping, split_files)
        sample_log['renames'] = renames
        
        for r in renames:
            icon = '✓' if r['status'] == 'SUCCESS' else '✗'
            print(f"    {icon} {r['original']} → {r['new']}")
        
        # ================================================================
        # STEP 5: Determine status
        # ================================================================
        successful_renames = [r for r in renames if r['status'] == 'SUCCESS']
        
        if len(successful_renames) == len(mapping):
            sample_log['status'] = 'SUCCESS'
            processing_log['summary']['success'] += 1
            successful_accessions[key].append(accession)
        elif successful_renames:
            sample_log['status'] = 'PARTIAL'
            processing_log['summary']['partial'] += 1
        else:
            sample_log['status'] = 'FAILED'
            processing_log['summary']['failed'] += 1
        
        print(f"  Final status: {sample_log['status']}")
        
        # List final files
        print(f"  Final files:")
        for f in sorted(os.listdir(sample_dir)):
            if f.endswith('.fastq.gz'):
                size = os.path.getsize(os.path.join(sample_dir, f))
                print(f"    {f} ({size:,} bytes)")
        
        processing_log['samples'].append(sample_log)


# ============================================================================
# SAVE LOG
# ============================================================================

log_file = os.path.join(log_dir, f'download_log_{datetime.now().strftime("%Y%m%d_%H%M%S")}.json')
with open(log_file, 'w') as f:
    json.dump(processing_log, f, indent=2, default=str)


# ============================================================================
# CREATE SUCCESSFUL-ONLY METADATA DICTIONARY
# ============================================================================

print(f"\n{'='*70}")
print("Updating metadata with successful samples only...")
print(f"{'='*70}")

# Create a filtered copy of gse_dict with only successful samples
gse_dict_successful = {}

for key in gse_dict.keys():
    if key in successful_accessions and len(successful_accessions[key]) > 0:
        # Deep copy the dataframe
        df_copy = gse_dict[key].copy()
        
        # Filter to only successful accessions
        df_filtered = df_copy[df_copy['run_accession'].isin(successful_accessions[key])]
        
        if len(df_filtered) > 0:
            gse_dict_successful[key] = df_filtered
            print(f"  {key}: {len(df_filtered)}/{len(gse_dict[key])} samples successful")
        else:
            print(f"  {key}: No successful samples (0/{len(gse_dict[key])})")
    else:
        print(f"  {key}: No successful samples (0/{len(gse_dict[key])})")

# Save successful-only dictionary
successful_dict_file = os.path.join(metadata_dir, 'dictionary_file_successful.pkl')
with open(successful_dict_file, 'wb') as f:
    pickle.dump(gse_dict_successful, f)
print(f"\nSaved: {successful_dict_file}")

# Also save as JSON for easy inspection
successful_json_file = os.path.join(metadata_dir, 'successful_samples.json')
successful_summary = {
    'timestamp': datetime.now().isoformat(),
    'bioprojects': {}
}
for key, df in gse_dict_successful.items():
    successful_summary['bioprojects'][key] = {
        'total_original': len(gse_dict[key]),
        'total_successful': len(df),
        'accessions': df['run_accession'].tolist()
    }

with open(successful_json_file, 'w') as f:
    json.dump(successful_summary, f, indent=2)
print(f"Saved: {successful_json_file}")


# ============================================================================
# SUMMARY
# ============================================================================

print(f"\n{'='*70}")
print("DOWNLOAD SUMMARY")
print(f"{'='*70}")

s = processing_log['summary']
print(f"\nTotal: {processing_log['total_samples']}")
print(f"  Success: {s['success']} (includes {s['skipped_existing']} skipped/existing)")
print(f"  Partial: {s['partial']}")
print(f"  Failed: {s['failed']}")

print(f"\nMapping method:")
print(f"  Positional (Entrez order): {s['used_positional_mapping']}")
print(f"  Content analysis (fallback): {s['used_content_analysis']}")

print(f"\nValidation:")
print(f"  Samples with warnings: {s['validation_warnings']}")

# Show samples with validation warnings
samples_with_warnings = [s for s in processing_log['samples'] if s.get('validation_warnings')]
if samples_with_warnings:
    print(f"\nSamples with validation warnings ({len(samples_with_warnings)}):")
    for s in samples_with_warnings[:5]:
        print(f"  {s['accession']}:")
        for w in s['validation_warnings']:
            print(f"    - {w}")
    if len(samples_with_warnings) > 5:
        print(f"  ... and {len(samples_with_warnings) - 5} more")

# Show failed samples
failed_samples = [s for s in processing_log['samples'] if s['status'] in ['DOWNLOAD_FAILED', 'FAILED']]
if failed_samples:
    print(f"\nFailed samples ({len(failed_samples)}):")
    for s in failed_samples[:10]:
        print(f"  {s['accession']}: {s['issues']}")
    if len(failed_samples) > 10:
        print(f"  ... and {len(failed_samples) - 10} more")

print(f"\nLog saved: {log_file}")
print(f"Successful metadata: {successful_dict_file}")


# ============================================================================
# QC - FastQC and MultiQC
# ============================================================================

print(f"\n{'='*70}")
print("Running FastQC...")
print(f"{'='*70}")

# Collect all FASTQ files for QC
qc_files = []
for key in gse_dict.keys():
    for accession in gse_dict[key]['run_accession']:
        sample_dir = os.path.join(output_dir, 'fastq', key, accession)
        if not os.path.exists(sample_dir):
            continue
        for f in os.listdir(sample_dir):
            if f.endswith(('.fastq.gz', '.fq.gz')):
                filepath = os.path.join(sample_dir, f)
                # Only include non-empty files
                if os.path.getsize(filepath) > 0:
                    qc_files.append(filepath)

print(f"Found {len(qc_files)} FASTQ files for QC")

# Run FastQC with multithreading
# FastQC -t flag specifies threads PER FILE, so we batch to control parallelism
if qc_files:
    # FastQC can process multiple files in parallel
    # -t specifies threads to use (each file gets 1 thread by default)
    # We'll process in batches, with each batch using available threads
    
    # Calculate optimal batch size based on available threads
    # FastQC uses ~1 thread per file when -t is specified
    batch_size = min(computing_threads, len(qc_files))
    
    total_batches = (len(qc_files) + batch_size - 1) // batch_size
    
    for i in range(0, len(qc_files), batch_size):
        batch = qc_files[i:i+batch_size]
        batch_num = (i // batch_size) + 1
        print(f"  FastQC batch {batch_num}/{total_batches}: {len(batch)} files...")
        
        # -t threads: number of files to process simultaneously
        fastqc_result = subprocess.run(
            ["fastqc", 
             "-t", str(min(computing_threads, len(batch))),  # Threads = number of parallel files
             "-o", qc_dir,
             "--quiet"] + batch,
            capture_output=True,
            text=True
        )
        
        if fastqc_result.returncode != 0:
            print(f"    Warning: FastQC returned non-zero exit code")
            if fastqc_result.stderr:
                print(f"    stderr: {fastqc_result.stderr[:200]}")
    
    print(f"FastQC complete. Reports in: {qc_dir}")
else:
    print("No FASTQ files found for QC!")

# ============================================================================
# Run MultiQC to aggregate reports
# ============================================================================

print(f"\n{'='*70}")
print("Running MultiQC...")
print(f"{'='*70}")

# MultiQC aggregates all FastQC reports into a single HTML report
# Note: MultiQC itself is single-threaded but is very fast
multiqc_report = os.path.join(qc_dir, 'multiqc_report.html')

# Remove old report if exists (to ensure fresh generation)
if os.path.exists(multiqc_report):
    os.remove(multiqc_report)
    print(f"  Removed old MultiQC report")

# Run MultiQC
print(f"  Aggregating FastQC reports...")
multiqc_result = subprocess.run(
    ["multiqc", 
     "--force",           # Overwrite existing reports
     "--no-data-dir",     # Don't create data directory (saves space)
     "-o", qc_dir,        # Output directory
     qc_dir],             # Input directory (where FastQC reports are)
    capture_output=True,
    text=True
)

if multiqc_result.returncode != 0:
    print(f"  Warning: MultiQC returned non-zero exit code: {multiqc_result.returncode}")
    if multiqc_result.stderr:
        print(f"  stderr: {multiqc_result.stderr[:500]}")
    if multiqc_result.stdout:
        print(f"  stdout: {multiqc_result.stdout[:500]}")

# Verify MultiQC report was created (CRITICAL for Snakemake)
if os.path.exists(multiqc_report):
    report_size = os.path.getsize(multiqc_report)
    print(f"  ✓ MultiQC report created: {multiqc_report} ({report_size:,} bytes)")
else:
    # If MultiQC failed, create a minimal placeholder report
    # This ensures Snakemake doesn't fail, but alerts user to the issue
    print(f"  ✗ MultiQC report not found! Creating placeholder...")
    
    placeholder_html = f"""<!DOCTYPE html>
<html>
<head><title>MultiQC Report - Generation Failed</title></head>
<body>
<h1>MultiQC Report Generation Failed</h1>
<p>Timestamp: {datetime.now().isoformat()}</p>
<p>The MultiQC report could not be generated automatically.</p>
<h2>Summary</h2>
<ul>
<li>Total samples processed: {processing_log['total_samples']}</li>
<li>Successful: {processing_log['summary']['success']}</li>
<li>Failed: {processing_log['summary']['failed']}</li>
<li>FASTQ files found: {len(qc_files)}</li>
</ul>
<h2>To generate manually:</h2>
<pre>cd {qc_dir}
multiqc .</pre>
<h2>Successful samples by BioProject:</h2>
<ul>
"""
    for key, accessions in successful_accessions.items():
        placeholder_html += f"<li>{key}: {len(accessions)} samples</li>\n"
    
    placeholder_html += """</ul>
</body>
</html>"""
    
    with open(multiqc_report, 'w') as f:
        f.write(placeholder_html)
    
    print(f"  Created placeholder report: {multiqc_report}")

# ============================================================================
# FINAL SUMMARY
# ============================================================================

print(f"\n{'='*70}")
print("PIPELINE COMPLETE")
print(f"{'='*70}")

print(f"\nOutputs:")
print(f"  FASTQ files: {output_dir}/fastq/")
print(f"  QC reports:  {qc_dir}/")
print(f"  MultiQC:     {multiqc_report}")
print(f"  Logs:        {log_dir}/")
print(f"  Metadata:    {metadata_dir}/")

print(f"\nKey files:")
print(f"  - dictionary_file_successful.pkl  (successful samples only)")
print(f"  - successful_samples.json         (human-readable summary)")
print(f"  - download_log_*.json             (detailed processing log)")

print(f"\n{'='*70}")
print("Done!")
print(f"{'='*70}")

sys.exit(0)
