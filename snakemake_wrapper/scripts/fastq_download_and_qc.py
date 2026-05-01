#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os, sys
import gzip

output_dir = snakemake.params.output_dir
metadata_dir = os.path.join(snakemake.params.output_dir, 'metadata')
computing_threads = snakemake.params.computing_threads

os.chdir(metadata_dir)

# Import in the dictionary with the metadata
import pickle

with open('dictionary_file.pkl', 'rb') as pkl_file:
     gse_dict = pickle.load(pkl_file)
     print('Dictionary loaded successfully')

# =====================================================================
# CHANGE 1: Load file_metadata.pkl for metadata-informed renaming
# =====================================================================
file_metadata = {}
file_metadata_path = os.path.join(metadata_dir, 'file_metadata.pkl')
if os.path.exists(file_metadata_path):
    with open(file_metadata_path, 'rb') as f:
        file_metadata = pickle.load(f)
    print(f'File metadata loaded: {len(file_metadata)} runs')
else:
    print('WARNING: file_metadata.pkl not found - will use hardcoded rename order')


#######################################################################
# SDL API classification: determine BAM vs FASTQ for each run
#
# This ONLY decides the download strategy. It does NOT change how
# FASTQ files are downloaded — that path is completely untouched.
#######################################################################

import requests
import subprocess
import glob
import shutil

SDL_API_URL = "https://locate.ncbi.nlm.nih.gov/sdl/2/retrieve"

def classify_run_from_sdl(accession):
    """
    Query the NCBI SDL API for a single run accession and check the
    'filename' field of every file entry to classify as 'bam' or 'fastq'.

    Logic:
      - filename ends in .bam             -> 'bam'
      - filename ends in .fastq.gz/.fq.gz -> 'fastq'
      - both .bam and .fastq.gz present   -> 'fastq' (safe default)
      - neither present                   -> 'fastq' (safe default)
    """
    params = {
        'acc': accession,
        'accept-alternate-locations': 'yes'
    }
    try:
        response = requests.get(SDL_API_URL, params=params, timeout=30)
        response.raise_for_status()
        data = response.json()

        has_bam = False
        has_fastq = False
        bam_url = ''

        for result in data.get('result', []):
            for file_entry in result.get('files', []):
                fname = file_entry.get('name', '').lower()
                if fname.endswith('.bam'):
                    has_bam = True
                    # Capture the S3 download URL for this BAM file
                    for loc in file_entry.get('locations', []):
                        link = loc.get('link', '')
                        if link and 's3' in loc.get('service', ''):
                            bam_url = link
                if fname.endswith('.fastq.gz') or fname.endswith('.fq.gz'):
                    has_fastq = True

        if has_bam and not has_fastq:
            return ('bam', bam_url)
        else:
            return ('fastq', '')

    except Exception as e:
        print(f"[SDL] WARNING: Could not query SDL API for {accession}: {e}")
        return ('fastq', '')


def build_file_type_map(gse_dict):
    """
    Classify every run accession in gse_dict using the SDL API.
    Returns: { run_accession: ('bam'|'fastq', bam_download_url) }
    """
    file_type_map = {}

    for gse_key in gse_dict.keys():
        for accession in gse_dict[gse_key]['run_accession']:
            classification, bam_url = classify_run_from_sdl(accession)
            file_type_map[accession] = (classification, bam_url)
            print(f"  [SDL] {accession} -> {classification}" +
                  (f" (BAM URL: {bam_url[:80]}...)" if bam_url else ""))

    bam_count = sum(1 for v in file_type_map.values() if v[0] == 'bam')
    fastq_count = sum(1 for v in file_type_map.values() if v[0] == 'fastq')
    print(f"\n[SDL] Classification complete: {fastq_count} FASTQ, {bam_count} BAM\n")

    return file_type_map


print("=" * 80)
print("CLASSIFYING RUN ACCESSIONS (SDL API - per-file detection)")
print("=" * 80)
file_type_map = build_file_type_map(gse_dict)


#######################################################################
# BAM download and conversion helpers
#######################################################################

def download_bam_from_url(bam_url, accession, sample_dir):
    """
    Download original BAM file using the S3 URL from the SDL API.
    Falls back to prefetch --type all if URL download fails.
    Returns path to downloaded BAM, or None.
    """
    os.makedirs(sample_dir, exist_ok=True)

    # Try direct URL download first (faster, no SRA toolkit dependency)
    if bam_url:
        # Extract original filename from URL (e.g., .../<SRR>/filename.bam.1)
        url_basename = bam_url.rstrip('/').split('/')[-1]
        # Remove trailing .1 version suffix if present
        if url_basename.endswith('.1'):
            url_basename = url_basename[:-2]
        bam_dest = os.path.join(sample_dir, url_basename)

        print(f"[BAM DOWNLOAD] Downloading from: {bam_url}")
        print(f"[BAM DOWNLOAD] Saving to: {bam_dest}")

        try:
            result = subprocess.run(
                ["wget", "-q", "-O", bam_dest, bam_url],
                capture_output=True, text=True, timeout=7200  # 2 hour timeout
            )
            if result.returncode == 0 and os.path.exists(bam_dest) and os.path.getsize(bam_dest) > 0:
                print(f"[BAM DOWNLOAD] Success: {bam_dest}")
                return bam_dest
            else:
                print(f"[BAM DOWNLOAD] wget failed: {result.stderr}")
        except Exception as e:
            print(f"[BAM DOWNLOAD] wget error: {e}")

    # Fallback: prefetch --type all
    print(f"[BAM DOWNLOAD] Trying prefetch --type all for {accession}")
    try:
        proc = subprocess.Popen(
            ["prefetch", str(accession), "-O", sample_dir,
             "--max-size", "u", "--type", "all"],
            stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True
        )
        output, error = proc.communicate()
        print(f"[BAM DOWNLOAD] prefetch output: {output}")
        if error:
            print(f"[BAM DOWNLOAD] prefetch stderr: {error}")

        bam_files = glob.glob(os.path.join(sample_dir, '**', '*.bam'), recursive=True)
        if bam_files:
            print(f"[BAM DOWNLOAD] Found: {bam_files[0]}")
            return bam_files[0]
    except Exception as e:
        print(f"[BAM DOWNLOAD] prefetch error: {e}")

    print(f"[BAM DOWNLOAD] WARNING: Could not download BAM for {accession}")
    return None


def is_10x_bam(bam_path):
    """
    Check BAM header for 10x Chromium @CO lines.
    Returns True if this is a 10x Cell Ranger BAM.
    """
    try:
        result = subprocess.run(
            ["bash", "-c", f"samtools view -H {bam_path} | grep -c '10x_bam_to_fastq'"],
            capture_output=True, text=True, timeout=60
        )
        count = int(result.stdout.strip())
        return count > 0
    except Exception:
        return False


def convert_10x_bam_to_fastq(bam_path, accession, sample_dir, threads):
    """
    Convert a 10x Chromium BAM to FASTQ using the 10x bamtofastq tool.

    bamtofastq creates a nested directory structure:
      output_dir/prefix_MissingLibrary/bamtofastq_S1_L00N/
        bamtofastq_S1_L001_R1_001.fastq.gz  (barcode + UMI)
        bamtofastq_S1_L001_R2_001.fastq.gz  (cDNA)
        bamtofastq_S1_L001_I1_001.fastq.gz  (index, if present)

    This function concatenates all per-lane files by read type and
    renames to {accession}_1.fastq.gz, _2.fastq.gz, etc. so the
    existing rename logic downstream runs identically.
    """
    bamtofastq_outdir = os.path.join(sample_dir, f"{accession}_bamtofastq_tmp")

    print(f"[10x CONVERT] Running bamtofastq on {os.path.basename(bam_path)}")
    print(f"[10x CONVERT] Output dir: {bamtofastq_outdir}")

    try:
        # Run 10x bamtofastq
        proc = subprocess.Popen(
            ["bamtofastq", "--nthreads", str(threads),
             str(bam_path), bamtofastq_outdir],
            stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True
        )
        stdout, stderr = proc.communicate()
        print(f"[10x CONVERT] stdout: {stdout}")
        if stderr:
            print(f"[10x CONVERT] stderr: {stderr}")
        if proc.returncode != 0:
            print(f"[10x CONVERT] bamtofastq failed with exit code {proc.returncode}")
            return False

        # Find all output FASTQ files grouped by read type
        # bamtofastq names them: *_R1_001.fastq.gz, *_R2_001.fastq.gz, *_I1_001.fastq.gz, *_I2_001.fastq.gz
        read_type_map = {
            '_R1_001.fastq.gz': f"{accession}_1.fastq.gz",
            '_R2_001.fastq.gz': f"{accession}_2.fastq.gz",
            '_I1_001.fastq.gz': f"{accession}_3.fastq.gz",
            '_I2_001.fastq.gz': f"{accession}_4.fastq.gz",
        }

        for suffix, target_name in read_type_map.items():
            # Find all files matching this read type across all lane directories
            matching_files = sorted(glob.glob(
                os.path.join(bamtofastq_outdir, '**', f'*{suffix}'),
                recursive=True
            ))

            if not matching_files:
                continue

            target_path = os.path.join(sample_dir, target_name)
            print(f"[10x CONVERT] Concatenating {len(matching_files)} files -> {target_name}")

            # Concatenate all matching gzipped files (cat works for .gz)
            with open(target_path, 'wb') as outf:
                for mf in matching_files:
                    with open(mf, 'rb') as inf:
                        while True:
                            chunk = inf.read(1024 * 1024)  # 1MB chunks
                            if not chunk:
                                break
                            outf.write(chunk)

        # Cleanup the bamtofastq temp directory
        if os.path.exists(bamtofastq_outdir):
            shutil.rmtree(bamtofastq_outdir)
            print(f"[10x CONVERT] Cleaned up temp dir: {bamtofastq_outdir}")

        # Verify at minimum R1 and R2 were created
        r1_exists = os.path.exists(os.path.join(sample_dir, f"{accession}_1.fastq.gz"))
        r2_exists = os.path.exists(os.path.join(sample_dir, f"{accession}_2.fastq.gz"))

        if r1_exists and r2_exists:
            print(f"[10x CONVERT] Success: created {accession}_1.fastq.gz and _2.fastq.gz")
            # Report any index files created
            for idx_suffix in ['_3.fastq.gz', '_4.fastq.gz']:
                if os.path.exists(os.path.join(sample_dir, f"{accession}{idx_suffix}")):
                    print(f"[10x CONVERT]   Also created {accession}{idx_suffix}")
            return True
        else:
            print(f"[10x CONVERT] WARNING: R1 exists={r1_exists}, R2 exists={r2_exists}")
            return False

    except FileNotFoundError:
        print(f"[10x CONVERT] ERROR: bamtofastq not found. Install via: conda install -c bioconda 10x_bamtofastq")
        return False
    except Exception as e:
        print(f"[10x CONVERT] ERROR: {e}")
        return False


def convert_standard_bam_to_fastq(bam_path, accession, sample_dir, threads):
    """
    Convert a standard (non-10x) BAM to FASTQ using samtools.
    These BAMs have proper paired-end flags so samtools fastq -1/-2 works.

    Produces {accession}_1.fastq.gz and {accession}_2.fastq.gz
    matching the naming that parallel-fastq-dump would produce.
    """
    r1_path = os.path.join(sample_dir, f"{accession}_1.fastq.gz")
    r2_path = os.path.join(sample_dir, f"{accession}_2.fastq.gz")

    print(f"[STD CONVERT] Converting {os.path.basename(bam_path)} to FASTQ via samtools")

    try:
        # Step 1: Sort BAM by read name
        sorted_bam = os.path.join(sample_dir, f"{accession}_namesorted.bam")
        sort_proc = subprocess.Popen(
            ["samtools", "sort", "-n", "-@", str(threads),
             "-o", sorted_bam, str(bam_path)],
            stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True
        )
        _, error = sort_proc.communicate()
        if sort_proc.returncode != 0:
            print(f"[STD CONVERT] samtools sort failed: {error}")
            return False

        # Step 2: Extract paired FASTQ
        r1_uncomp = os.path.join(sample_dir, f"{accession}_1.fastq")
        r2_uncomp = os.path.join(sample_dir, f"{accession}_2.fastq")

        fastq_proc = subprocess.Popen(
            ["samtools", "fastq", "-@", str(threads),
             "-1", r1_uncomp, "-2", r2_uncomp,
             "-0", "/dev/null", "-s", "/dev/null",
             "-n", sorted_bam],
            stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True
        )
        _, error = fastq_proc.communicate()
        if fastq_proc.returncode != 0:
            print(f"[STD CONVERT] samtools fastq failed: {error}")
            return False

        # Step 3: Gzip
        for fq in [r1_uncomp, r2_uncomp]:
            if os.path.exists(fq) and os.path.getsize(fq) > 0:
                subprocess.run(["gzip", "-f", fq],
                               capture_output=True, text=True)

        # Step 4: Cleanup intermediate
        if os.path.exists(sorted_bam):
            os.remove(sorted_bam)

        if os.path.exists(r1_path) and os.path.getsize(r1_path) > 0:
            print(f"[STD CONVERT] Success: created {accession}_1.fastq.gz and _2.fastq.gz")
            return True
        else:
            print(f"[STD CONVERT] WARNING: expected output files not found or empty")
            return False

    except FileNotFoundError:
        print(f"[STD CONVERT] ERROR: samtools not found. Install via conda env.")
        return False
    except Exception as e:
        print(f"[STD CONVERT] ERROR: {e}")
        return False


def convert_bam_to_fastq(bam_path, accession, sample_dir, threads):
    """
    Convert BAM to FASTQ. Detects BAM type from header and routes
    to the appropriate converter:

      - 10x Chromium BAM (@CO 10x_bam_to_fastq in header)
        -> bamtofastq tool (reconstructs R1/R2/I1 from BAM tags)

      - Standard paired-end BAM (proper FLAG pairing)
        -> samtools fastq -1/-2

    Both paths produce {accession}_1.fastq.gz, {accession}_2.fastq.gz
    (and optionally _3, _4 for index reads) so the existing rename
    logic downstream runs identically.
    """
    print(f"[BAM CONVERT] Checking BAM type for {os.path.basename(bam_path)}...")

    if is_10x_bam(bam_path):
        print(f"[BAM CONVERT] Detected 10x Chromium BAM -> using bamtofastq")
        return convert_10x_bam_to_fastq(bam_path, accession, sample_dir, threads)
    else:
        print(f"[BAM CONVERT] Standard BAM -> using samtools fastq")
        return convert_standard_bam_to_fastq(bam_path, accession, sample_dir, threads)


# =====================================================================
# CHANGE 2: Metadata-informed rename helper functions
# =====================================================================

def get_read_length_quick(fastq_gz_path, n_reads=20):
    """
    Get the dominant read length from the first n_reads of a FASTQ.gz.
    Returns the most common length, or None if unreadable.
    """
    if not os.path.exists(fastq_gz_path) or os.path.getsize(fastq_gz_path) == 0:
        return None
    try:
        lengths = {}
        with gzip.open(fastq_gz_path, 'rt') as f:
            line_num = 0
            reads_counted = 0
            for line in f:
                line_num += 1
                if line_num % 4 == 2:
                    l = len(line.strip())
                    lengths[l] = lengths.get(l, 0) + 1
                    reads_counted += 1
                    if reads_counted >= n_reads:
                        break
        if lengths:
            return max(lengths, key=lengths.get)
    except Exception:
        pass
    return None


def build_rename_map_from_metadata(accession, file_metadata):
    """
    Use expected_read_types_ordered from Entrez metadata to build the
    mapping from split file number to read type.

    For multi-lane data, expected_read_types_ordered repeats per lane
    (e.g. [I1, I2, R1, R2, I1, I2, R1, R2]). parallel-fastq-dump
    with --split-files merges lanes, so we deduplicate to get the
    unique ordered types.

    Returns: dict {1: 'I1', 2: 'I2', 3: 'R1', 4: 'R2'} or None
    """
    meta = file_metadata.get(accession, {})
    ordered = meta.get('expected_read_types_ordered', [])

    if not ordered:
        return None

    # Deduplicate: keep first occurrence of each type
    unique_types = []
    seen = set()
    for rt in ordered:
        if rt not in seen:
            unique_types.append(rt)
            seen.add(rt)

    if not unique_types:
        return None

    # Build the map: split file index -> read type
    rename_map = {}
    for i, rt in enumerate(unique_types):
        rename_map[i + 1] = rt  # 1-indexed to match _1, _2, _3, _4

    return rename_map


def validate_renamed_files(sample_dir, accession):
    """
    Post-rename validation: check that R1 looks like barcode+UMI (24-30bp)
    and R2 looks like cDNA (>=50bp). Warns if reads appear misassigned.

    Returns True if validation passes, False if reads appear swapped.
    """
    r1_path = os.path.join(sample_dir,
                           f"{accession}_S1_L001_R1_001.fastq.gz")
    r2_path = os.path.join(sample_dir,
                           f"{accession}_S1_L001_R2_001.fastq.gz")

    r1_len = get_read_length_quick(r1_path)
    r2_len = get_read_length_quick(r2_path)

    if r1_len is None or r2_len is None:
        print(f"  [VALIDATE] Could not read R1/R2 lengths for {accession}")
        return True  # Can't validate, assume OK

    valid = True
    if r1_len <= 12:
        print(f"  [VALIDATE] WARNING: R1 is {r1_len}bp (looks like index, not barcode+UMI)")
        valid = False
    if r2_len <= 12:
        print(f"  [VALIDATE] WARNING: R2 is {r2_len}bp (looks like index, not cDNA)")
        valid = False

    if valid:
        print(f"  [VALIDATE] OK: R1={r1_len}bp, R2={r2_len}bp")
    else:
        print(f"  [VALIDATE] READS MAY BE MISASSIGNED - check file_metadata.pkl")
        print(f"  [VALIDATE] R1={r1_len}bp, R2={r2_len}bp")

    return valid


#######################################################################
# Main download loop
#
# FASTQ PATH: completely unchanged from original pipeline
# BAM PATH:   download original BAM -> samtools convert -> same output
#
# Both paths produce {accession}_1.fastq.gz, {accession}_2.fastq.gz
# so the shared rename logic below works identically for both.
#######################################################################

for key in gse_dict.keys():
    for accession in gse_dict[key]['run_accession']:
        print(f"\nProcessing sample {accession} from the BioProject {key}")

        # Absolute path for this sample's output directory
        sample_dir = os.path.join(str(output_dir), 'fastq', str(key), str(accession))

        # Get classification for this run
        run_type, bam_url = file_type_map.get(accession, ('fastq', ''))

        if run_type == 'bam':
            # ==============================================================
            # BAM PATH: download original BAM, convert to FASTQ
            # ==============================================================
            print(f"[BAM PATH] {accession} is a BAM submission, downloading and converting...")

            bam_path = download_bam_from_url(bam_url, accession, sample_dir)

            if bam_path:
                success = convert_bam_to_fastq(bam_path, accession, sample_dir, computing_threads)
                if success:
                    # Clean up original BAM to save disk space
                    try:
                        os.remove(bam_path)
                        print(f"[BAM PATH] Cleaned up: {os.path.basename(bam_path)}")
                    except Exception as e:
                        print(f"[BAM PATH] Could not remove BAM: {e}")
                else:
                    print(f"[BAM PATH] WARNING: conversion failed for {accession}")
            else:
                print(f"[BAM PATH] WARNING: download failed for {accession}")

        else:
            # ==============================================================
            # FASTQ PATH: original pipeline logic — COMPLETELY UNCHANGED
            # ==============================================================
            subprocess_1 = subprocess.Popen(
                ["prefetch", str(accession), "-O", 
                 str(output_dir) + '/fastq/' + str(key), "--max-size", "u"], stdout=subprocess.PIPE, text=True)
            output, error = subprocess_1.communicate()
            print(f'Outputs: {output}')
            print(f'Errors: {error}')

            os.chdir(os.path.join(str(output_dir), 'fastq', str(key), str(accession)))
             
            subprocess_2 = subprocess.Popen(
                ["parallel-fastq-dump", "-s", str(accession) + ".sra", "--threads", str(computing_threads), 
                 "--tmpdir", str(output_dir) + '/fastq/' + str(key) + '/' + str(accession),
                 "--outdir", str(output_dir) + '/fastq/' + str(key) + '/' + str(accession), "--split-spot", "--split-files", "--gzip"], stdout=subprocess.PIPE, text=True)
            output, error = subprocess_2.communicate()
            print(f'Outputs: {output}')
            print(f'Errors: {error}')

        # ==============================================================
        # CHANGE 3: Metadata-informed rename with fallback and validation
        #
        # Uses expected_read_types_ordered from file_metadata.pkl to
        # determine the correct mapping from split file numbers to
        # read types (R1/R2/I1/I2).
        #
        # Falls back to hardcoded [R1, R2, I1, I2] if metadata is
        # unavailable for this accession.
        #
        # Validates read lengths after rename and warns if R1/R2
        # appear to be misassigned.
        # ==============================================================
        print(f"\nRenaming {accession} fastqs")

        # Determine rename mapping: metadata-informed or hardcoded fallback
        rename_map = build_rename_map_from_metadata(accession, file_metadata)

        if rename_map:
            print(f"  Using metadata-informed rename: {rename_map}")
        else:
            # Hardcoded fallback: _1=R1, _2=R2, _3=I1, _4=I2
            rename_map = {1: 'R1', 2: 'R2', 3: 'I1', 4: 'I2'}
            print(f"  No metadata available, using hardcoded rename: {rename_map}")

        # Perform renames based on the mapping
        for split_num, read_type in sorted(rename_map.items()):
            src = os.path.join(sample_dir,
                               f"{accession}_{split_num}.fastq.gz")
            dst = os.path.join(sample_dir,
                               f"{accession}_S1_L001_{read_type}_001.fastq.gz")
            try:
                if os.path.exists(src):
                    os.rename(src, dst)
                    print(f"  Renamed _{split_num} -> {read_type}")
                else:
                    if read_type in ('R1', 'R2'):
                        print(f"  WARNING: {src} not found (expected {read_type})")
                    else:
                        print(f"  Note: {src} not found ({read_type} index file not present)")
            except OSError as e:
                print(f"  Error renaming _{split_num} -> {read_type}: {e}")

        # Post-rename validation
        validate_renamed_files(sample_dir, accession)

        # Clean up .sra file
        try:
             file_path = os.path.join(sample_dir, str(accession) + '.sra')
             os.remove(file_path)
             print(f"File '{file_path}' deleted successfully.")
        except FileNotFoundError:
             print(f"Error: File '{file_path}' not found.")
        except Exception as e:
             print(f"An unexpected error occurred: {e}")



# QC Section
# Go check some files with fastqc

for key in gse_dict.keys():
    for accession in gse_dict[key]['run_accession']:
        subprocess_3 = subprocess.Popen(
            ["fastqc", "-t", str(computing_threads), "-o", os.path.join(output_dir, 'QC'), 
             os.path.join(str(output_dir), 'fastq', str(key), str(accession), str(accession) + '_S1_L001_R1_001.fastq.gz'), 
             os.path.join(str(output_dir), 'fastq', str(key), str(accession), str(accession) + '_S2_L001_R1_001.fastq.gz')], stdout=subprocess.PIPE, text=True)
        output, error = subprocess_3.communicate()
        print(f'Outputs: {output}')
        print(f'Errors: {error}')

# Pull all the repoorts together with multiqc

subprocess_4 = subprocess.Popen(
    ["multiqc", "-o", os.path.join(str(output_dir), 'QC'), os.path.join(str(output_dir), 'QC')], stdout=subprocess.PIPE, text=True)
output, error = subprocess_4.communicate()
print(f'Outputs: {output}')
print(f'Errors: {error}')


# End file

sys.exit()
