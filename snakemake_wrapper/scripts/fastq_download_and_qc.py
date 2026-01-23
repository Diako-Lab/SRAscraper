#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys
import gzip
import subprocess
import pickle
import json
import re
from collections import Counter
from datetime import datetime
from pathlib import Path

# Snakemake parameters
output_dir = snakemake.params.output_dir
metadata_dir = os.path.join(snakemake.params.output_dir, 'metadata')
computing_threads = snakemake.params.computing_threads

os.chdir(metadata_dir)


# ============================================================================
# CONFIGURATION
# ============================================================================

# Read length thresholds for 10x content classification
BARCODE_UMI_LENGTHS = {24, 26, 28}  # R1 in 10x
SAMPLE_INDEX_MAX_LENGTH = 12        # I1/I2 in 10x
CDNA_MIN_LENGTH = 50                # R2 in 10x
SAMPLE_SIZE = 500


# ============================================================================
# CONTENT ANALYSIS FUNCTIONS
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
    """Classify 10x read type by length"""
    if length is None:
        return 'UNKNOWN'
    if length in BARCODE_UMI_LENGTHS:
        return 'BARCODE_UMI'  # This is R1 in 10x
    if length <= SAMPLE_INDEX_MAX_LENGTH:
        return 'SAMPLE_INDEX'  # This is I1/I2 in 10x
    if length >= CDNA_MIN_LENGTH:
        return 'CDNA'  # This is R2 in 10x
    return 'UNKNOWN'


def map_content_to_read_type(content_type):
    """Map content classification to standard read type"""
    mapping = {
        'BARCODE_UMI': 'R1',
        'CDNA': 'R2',
        'SAMPLE_INDEX': 'I',  # Will be I1, I2, etc.
    }
    return mapping.get(content_type)


def analyze_split_files(sample_dir, accession):
    """
    Analyze all split files from fastq-dump.
    Returns dict mapping filename to content analysis.
    """
    results = {}
    
    # Find split files (_1, _2, _3, _4, etc.)
    for f in os.listdir(sample_dir):
        if f.startswith(accession) and f.endswith('.fastq.gz'):
            # Check if it's a numbered split file
            match = re.match(rf'{accession}_(\d+)\.fastq\.gz', f)
            if match:
                split_num = int(match.group(1))
                filepath = os.path.join(sample_dir, f)
                
                length_info = sample_fastq_lengths(filepath)
                content_type = classify_read_by_length(length_info.get('dominant_length'))
                
                results[f] = {
                    'split_number': split_num,
                    'filepath': filepath,
                    'dominant_length': length_info.get('dominant_length'),
                    'content_type': content_type,
                    'read_type': map_content_to_read_type(content_type)
                }
    
    return results


def rename_by_content_analysis(sample_dir, accession, file_analysis, expected_types=None):
    """
    Rename split files based on content analysis.
    
    If expected_types is provided (from metadata), validate against it.
    
    Returns list of rename operations.
    """
    renames = []
    
    # Group files by content type
    content_groups = {
        'BARCODE_UMI': [],
        'CDNA': [],
        'SAMPLE_INDEX': [],
        'UNKNOWN': []
    }
    
    for filename, info in file_analysis.items():
        content_type = info.get('content_type', 'UNKNOWN')
        content_groups[content_type].append((filename, info))
    
    # Assign read types
    # R1 = BARCODE_UMI (24-28bp)
    for filename, info in content_groups['BARCODE_UMI']:
        new_name = f"{accession}_S1_L001_R1_001.fastq.gz"
        renames.append(create_rename_record(sample_dir, filename, new_name, 'R1', info))
    
    # R2 = CDNA (50+ bp)
    for filename, info in content_groups['CDNA']:
        new_name = f"{accession}_S1_L001_R2_001.fastq.gz"
        renames.append(create_rename_record(sample_dir, filename, new_name, 'R2', info))
    
    # I1, I2 = SAMPLE_INDEX (<=12bp)
    for i, (filename, info) in enumerate(content_groups['SAMPLE_INDEX']):
        idx_num = i + 1
        new_name = f"{accession}_S1_L001_I{idx_num}_001.fastq.gz"
        renames.append(create_rename_record(sample_dir, filename, new_name, f'I{idx_num}', info))
    
    # Perform renames
    for rename in renames:
        old_path = rename['old_path']
        new_path = rename['new_path']
        
        try:
            if old_path != new_path and os.path.exists(old_path):
                os.rename(old_path, new_path)
                rename['status'] = 'SUCCESS'
            elif old_path == new_path:
                rename['status'] = 'NO_CHANGE'
            else:
                rename['status'] = 'FILE_NOT_FOUND'
        except Exception as e:
            rename['status'] = 'ERROR'
            rename['error'] = str(e)
    
    return renames


def create_rename_record(sample_dir, old_name, new_name, read_type, info):
    """Create a rename record dict"""
    return {
        'original': old_name,
        'new': new_name,
        'read_type': read_type,
        'content_type': info.get('content_type'),
        'length': info.get('dominant_length'),
        'old_path': os.path.join(sample_dir, old_name),
        'new_path': os.path.join(sample_dir, new_name),
        'status': None
    }


def rename_by_metadata_patterns(sample_dir, accession, file_analysis, expected_types):
    """
    Rename files when metadata tells us exactly what read types exist.
    
    Uses content analysis to MAP split files to expected types.
    Example: metadata says [I1, R1, R2], content says _1=28bp, _2=10bp, _3=90bp
    -> 28bp matches R1, 10bp matches I1, 90bp matches R2
    """
    renames = []
    
    # Build content-to-expected mapping
    # For 10x: BARCODE_UMI->R1, CDNA->R2, SAMPLE_INDEX->I1/I2
    content_to_expected = {}
    
    if 'R1' in expected_types:
        content_to_expected['BARCODE_UMI'] = 'R1'
    if 'R2' in expected_types:
        content_to_expected['CDNA'] = 'R2'
    
    # Count expected sample indices
    index_count = sum(1 for t in expected_types if t.startswith('I'))
    index_assigned = 0
    
    for filename, info in file_analysis.items():
        content_type = info.get('content_type')
        
        if content_type in content_to_expected:
            read_type = content_to_expected[content_type]
        elif content_type == 'SAMPLE_INDEX' and index_count > 0:
            index_assigned += 1
            read_type = f'I{index_assigned}'
        else:
            read_type = None
        
        if read_type:
            new_name = f"{accession}_S1_L001_{read_type}_001.fastq.gz"
            renames.append(create_rename_record(sample_dir, filename, new_name, read_type, info))
    
    # Perform renames
    for rename in renames:
        old_path = rename['old_path']
        new_path = rename['new_path']
        
        try:
            if old_path != new_path and os.path.exists(old_path):
                os.rename(old_path, new_path)
                rename['status'] = 'SUCCESS'
            elif old_path == new_path:
                rename['status'] = 'NO_CHANGE'
            else:
                rename['status'] = 'FILE_NOT_FOUND'
        except Exception as e:
            rename['status'] = 'ERROR'
            rename['error'] = str(e)
    
    return renames


def download_bam_directly(accession, sample_dir, meta):
    """
    Download BAM file directly from cloud storage.
    Uses URLs from metadata if available.
    """
    result = {
        'method': 'direct_bam',
        'files': [],
        'status': 'FAILED'
    }
    
    # Try to get BAM URLs from original files
    for f in meta.get('original_files', []):
        if f.get('is_bam') and f.get('url'):
            url = f['url']
            filename = f['filename']
            output_path = os.path.join(sample_dir, filename)
            
            try:
                # Handle different URL schemes
                if url.startswith('s3://'):
                    # Use AWS CLI
                    cmd = ['aws', 's3', 'cp', '--no-sign-request', url, output_path]
                elif url.startswith('gs://'):
                    # Use gsutil
                    cmd = ['gsutil', 'cp', url, output_path]
                else:
                    # Use wget for FTP/HTTP
                    cmd = ['wget', '-q', '-O', output_path, url]
                
                download_result = subprocess.run(cmd, capture_output=True, text=True, timeout=7200)
                
                if download_result.returncode == 0 and os.path.exists(output_path):
                    result['files'].append({
                        'filename': filename,
                        'path': output_path,
                        'status': 'SUCCESS'
                    })
                    result['status'] = 'SUCCESS'
                else:
                    result['files'].append({
                        'filename': filename,
                        'status': 'FAILED',
                        'error': download_result.stderr
                    })
            except Exception as e:
                result['files'].append({
                    'filename': filename,
                    'status': 'ERROR',
                    'error': str(e)
                })
    
    return result


# ============================================================================
# MAIN DOWNLOAD LOGIC
# ============================================================================

# Load metadata
with open('dictionary_file.pkl', 'rb') as f:
    gse_dict = pickle.load(f)
    print('Loaded: dictionary_file.pkl')

file_metadata = {}
if os.path.exists('file_metadata.pkl'):
    with open('file_metadata.pkl', 'rb') as f:
        file_metadata = pickle.load(f)
        print(f'Loaded: file_metadata.pkl ({len(file_metadata)} runs)')
else:
    print('WARNING: No file_metadata.pkl - will use content analysis for all files')

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
        'used_metadata_patterns': 0,
        'used_content_analysis': 0,
        'bam_downloads': 0,
        'kept_original': 0
    }
}


for key in gse_dict.keys():
    print(f"\n{'='*70}")
    print(f"Processing BioProject: {key}")
    print(f"{'='*70}")
    
    for accession in gse_dict[key]['run_accession']:
        processing_log['total_samples'] += 1
        
        sample_log = {
            'accession': accession,
            'bioproject': key,
            'status': None,
            'method': None,
            'expected_types': [],
            'renames': [],
            'issues': []
        }
        
        print(f"\n--- {accession} ---")
        
        sample_dir = os.path.join(output_dir, 'fastq', key, accession)
        os.makedirs(sample_dir, exist_ok=True)
        
        # Get metadata for this run
        meta = file_metadata.get(accession, {})
        expected_types = meta.get('file_types_present', [])
        is_bam = meta.get('is_bam', False)
        is_10x = meta.get('is_10x', False)
        needs_content_analysis = meta.get('needs_content_analysis', True)
        
        sample_log['expected_types'] = expected_types
        
        print(f"  Metadata: types={expected_types}, is_10x={is_10x}, is_bam={is_bam}")
        
        # ================================================================
        # Handle BAM files
        # ================================================================
        if is_bam:
            print(f"  Downloading BAM directly...")
            bam_result = download_bam_directly(accession, sample_dir, meta)
            
            if bam_result['status'] == 'SUCCESS':
                sample_log['status'] = 'SUCCESS'
                sample_log['method'] = 'BAM_DIRECT'
                processing_log['summary']['bam_downloads'] += 1
                processing_log['summary']['success'] += 1
                print(f"  ✓ BAM downloaded")
                for f in bam_result['files']:
                    print(f"    - {f['filename']}")
            else:
                print(f"  BAM download failed, trying fastq-dump fallback...")
                is_bam = False  # Fall through to fastq-dump
        
        # ================================================================
        # Standard SRA download with fastq-dump
        # ================================================================
        if not is_bam:
            # Step 1: Prefetch
            print(f"  Prefetching...")
            prefetch_result = subprocess.run(
                ["prefetch", accession, "-O",
                 os.path.join(output_dir, 'fastq', key), "--max-size", "u"],
                capture_output=True, text=True
            )
            
            if prefetch_result.returncode != 0:
                print(f"  ERROR: Prefetch failed")
                sample_log['status'] = 'PREFETCH_FAILED'
                sample_log['issues'].append(prefetch_result.stderr)
                processing_log['samples'].append(sample_log)
                processing_log['summary']['failed'] += 1
                continue
            
            # Step 2: Extract with parallel-fastq-dump
            print(f"  Extracting FASTQs...")
            os.chdir(sample_dir)
            
            sra_file = os.path.join(sample_dir, f"{accession}.sra")
            if not os.path.exists(sra_file):
                alt_sra = os.path.join(sample_dir, accession, f"{accession}.sra")
                if os.path.exists(alt_sra):
                    sra_file = alt_sra
            
            if not os.path.exists(sra_file):
                print(f"  ERROR: SRA file not found")
                sample_log['status'] = 'SRA_NOT_FOUND'
                processing_log['samples'].append(sample_log)
                processing_log['summary']['failed'] += 1
                continue
            
            extract_result = subprocess.run(
                ["parallel-fastq-dump", "-s", sra_file,
                 "--threads", str(computing_threads),
                 "--tmpdir", sample_dir,
                 "--outdir", sample_dir,
                 "--split-spot", "--split-files", "--gzip"],
                capture_output=True, text=True
            )
            
            if extract_result.returncode != 0:
                print(f"  ERROR: Extraction failed")
                sample_log['status'] = 'EXTRACTION_FAILED'
                sample_log['issues'].append(extract_result.stderr)
                processing_log['samples'].append(sample_log)
                processing_log['summary']['failed'] += 1
                continue
            
            # Step 3: Analyze split files
            print(f"  Analyzing files...")
            file_analysis = analyze_split_files(sample_dir, accession)
            
            for fname, info in file_analysis.items():
                print(f"    {fname}: {info['dominant_length']}bp → {info['content_type']}")
            
            # Step 4: Decide renaming strategy
            if expected_types and not needs_content_analysis:
                # We know from metadata what read types should exist
                print(f"  Renaming by metadata patterns: {expected_types}")
                renames = rename_by_metadata_patterns(sample_dir, accession, file_analysis, expected_types)
                sample_log['method'] = 'METADATA_PATTERNS'
                processing_log['summary']['used_metadata_patterns'] += 1
                
            elif is_10x or any(info['content_type'] in ['BARCODE_UMI', 'CDNA'] 
                              for info in file_analysis.values()):
                # Looks like 10x data, use content analysis
                print(f"  Renaming by content analysis (10x detected)")
                renames = rename_by_content_analysis(sample_dir, accession, file_analysis)
                sample_log['method'] = 'CONTENT_ANALYSIS'
                processing_log['summary']['used_content_analysis'] += 1
                
            else:
                # Unknown/bulk data - keep original names
                print(f"  Keeping original names (no patterns detected)")
                renames = []
                sample_log['method'] = 'KEEP_ORIGINAL'
                sample_log['issues'].append("No R1/R2 patterns - keeping original names")
                processing_log['summary']['kept_original'] += 1
            
            sample_log['renames'] = renames
            
            # Print rename results
            for r in renames:
                icon = '✓' if r['status'] == 'SUCCESS' else '·'
                print(f"    {icon} {r['original']} → {r['new']} ({r.get('length')}bp)")
            
            # Step 5: Cleanup SRA file
            try:
                if os.path.exists(sra_file):
                    os.remove(sra_file)
            except:
                pass
            
            # Determine final status
            if sample_log['issues']:
                sample_log['status'] = 'PARTIAL'
                processing_log['summary']['partial'] += 1
            else:
                sample_log['status'] = 'SUCCESS'
                processing_log['summary']['success'] += 1
        
        processing_log['samples'].append(sample_log)
        print(f"  Final: {sample_log['status']} ({sample_log['method']})")


# ============================================================================
# SAVE LOG
# ============================================================================

log_file = os.path.join(log_dir, f'download_log_{datetime.now().strftime("%Y%m%d_%H%M%S")}.json')
with open(log_file, 'w') as f:
    json.dump(processing_log, f, indent=2, default=str)


# ============================================================================
# SUMMARY
# ============================================================================

print(f"\n{'='*70}")
print("DOWNLOAD SUMMARY")
print(f"{'='*70}")

s = processing_log['summary']
print(f"\nTotal: {processing_log['total_samples']}")
print(f"  Success: {s['success']}")
print(f"  Partial: {s['partial']}")
print(f"  Failed: {s['failed']}")

print(f"\nMethods used:")
print(f"  Metadata patterns (R1/R2 from original names): {s['used_metadata_patterns']}")
print(f"  Content analysis (read length): {s['used_content_analysis']}")
print(f"  BAM direct download: {s['bam_downloads']}")
print(f"  Kept original names: {s['kept_original']}")

# Show samples needing review
needs_review = [s for s in processing_log['samples'] if s['method'] == 'KEEP_ORIGINAL']
if needs_review:
    print(f"\nSamples kept with original names ({len(needs_review)}):")
    for s in needs_review[:5]:
        print(f"  {s['accession']}")
    if len(needs_review) > 5:
        print(f"  ... and {len(needs_review) - 5} more")

print(f"\nLog saved: {log_file}")


# ============================================================================
# QC
# ============================================================================

print(f"\n{'='*70}")
print("Running FastQC...")
print(f"{'='*70}")

qc_files = []
for key in gse_dict.keys():
    for accession in gse_dict[key]['run_accession']:
        sample_dir = os.path.join(output_dir, 'fastq', key, accession)
        if not os.path.exists(sample_dir):
            continue
        for f in os.listdir(sample_dir):
            if f.endswith(('.fastq.gz', '.fq.gz')):
                qc_files.append(os.path.join(sample_dir, f))

if qc_files:
    batch_size = 20
    for i in range(0, len(qc_files), batch_size):
        batch = qc_files[i:i+batch_size]
        subprocess.run(
            ["fastqc", "-t", str(computing_threads), "-o", qc_dir] + batch,
            capture_output=True
        )
    subprocess.run(["multiqc", "-o", qc_dir, qc_dir], capture_output=True)
    print(f"QC reports: {qc_dir}")

print(f"\n{'='*70}")
print("Download complete!")
print(f"{'='*70}")

sys.exit(0)
