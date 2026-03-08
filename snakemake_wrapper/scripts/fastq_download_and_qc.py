#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os, sys

output_dir = snakemake.params.output_dir
metadata_dir = os.path.join(snakemake.params.output_dir, 'metadata')
computing_threads = snakemake.params.computing_threads

os.chdir(metadata_dir)

# Import in the dictionary with the metadata
import pickle

with open('dictionary_file.pkl', 'rb') as pkl_file:
     gse_dict = pickle.load(pkl_file)
     print('Dictionary loaded successfully')


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
        # SHARED: Rename files to Illumina convention (both paths)
        # ==============================================================
        print(f"\nRenamming {accession} fastqs")
        try:
             os.rename(os.path.join(str(output_dir), 'fastq', str(key), str(accession), str(accession) + '_4.fastq.gz'),
                       os.path.join(str(output_dir), 'fastq', str(key), str(accession), str(accession) + '_S1_L001_I2_001.fastq.gz'))
             print(f"\nIndex {accession}_4.fastq.gz file found in GEO repo. Confirm dual indexing in project.")
        except OSError as e:
            print("Error renaming file:", e)
        try:
             os.rename(os.path.join(str(output_dir), 'fastq', str(key), str(accession), str(accession) + '_3.fastq.gz'),
                       os.path.join(str(output_dir), 'fastq', str(key), str(accession), str(accession) + '_S1_L001_I1_001.fastq.gz'))
        except FileNotFoundError:
            print(f"\nIndex {accession}_3.fastq.gz file not found in GEO repo. Continuing with normal R1 R2 renaming.")
        except OSError as e:
            print("Error renaming file:", e)
        try:
            os.rename(os.path.join(str(output_dir), 'fastq', str(key), str(accession), str(accession) + '_1.fastq.gz'),
                      os.path.join(str(output_dir), 'fastq', str(key), str(accession), str(accession) + '_S1_L001_R1_001.fastq.gz'))          
            os.rename(os.path.join(str(output_dir), 'fastq', str(key), str(accession), str(accession) + '_2.fastq.gz'),
                      os.path.join(str(output_dir), 'fastq', str(key), str(accession), str(accession) + '_S1_L001_R2_001.fastq.gz')) 
        except FileNotFoundError:
            print(f"\nIndex {accession}_1/2.fastq.gz files not found in GEO repo.")
        except OSError as e:
            print("Error renaming file:", e)
        try:
             file_path = os.path.join(str(output_dir), 'fastq', str(key), str(accession), str(accession) + '.sra')
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
