#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os, sys, glob, subprocess

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
# STEP 1: Determine file type (BAM vs FASTQ) for every run accession
#         by inspecting SRA metadata BEFORE any downloads begin.
#######################################################################

from pysradb.sraweb import SRAweb
import pandas as pd

def build_file_type_map(gse_dict):
    """
    For every GSE project in gse_dict, fetch detailed SRA metadata
    and inspect the public_url field to classify each run accession
    as either 'bam' or 'fastq'.

    Returns a dict: { run_accession: 'bam' | 'fastq' }
    """
    sradb = SRAweb()
    file_type_map = {}

    # Collect all unique GSE keys to look up their SRP projects
    for gse_key in gse_dict.keys():
        print(f"\n[FILE TYPE CHECK] Looking up SRA metadata for {gse_key}...")
        try:
            gse_to_srp = sradb.gse_to_srp(gse_key)
            srp_ids = gse_to_srp['study_accession'].tolist()
        except Exception as e:
            print(f"[FILE TYPE CHECK] Could not map {gse_key} to SRP: {e}")
            # Default all runs in this project to fastq (existing behavior)
            for acc in gse_dict[gse_key]['run_accession']:
                file_type_map[acc] = 'fastq'
            continue

        for srp_id in srp_ids:
            try:
                metadata_df = sradb.sra_metadata(srp_id, detailed=True)

                for _, row in metadata_df.iterrows():
                    run_acc = row.get('run_accession', None)
                    if run_acc is None:
                        continue

                    # Scan ALL url columns for both .bam and .fastq patterns.
                    # Only classify as BAM when BAM urls are found and NO
                    # FASTQ urls are found.  If both exist, default to FASTQ
                    # so the existing prefetch + parallel-fastq-dump path is
                    # used — that path already works and we don't want to
                    # break it.
                    has_bam = False
                    has_fastq = False
                    for url_col in ['public_url', 'aws_url', 'gcp_url', 'ncbi_url']:
                        url_val = str(row.get(url_col, '')).lower()
                        if '.bam' in url_val:
                            has_bam = True
                        if '.fastq' in url_val or '.fq' in url_val:
                            has_fastq = True

                    if has_bam and not has_fastq:
                        file_type_map[run_acc] = 'bam'
                        print(f"  {run_acc} -> BAM (only BAM urls found, no FASTQ urls)")
                    elif has_bam and has_fastq:
                        file_type_map[run_acc] = 'fastq'
                        print(f"  {run_acc} -> FASTQ (both BAM and FASTQ urls found, defaulting to FASTQ)")
                    else:
                        file_type_map[run_acc] = 'fastq'
                        print(f"  {run_acc} -> FASTQ")

            except Exception as e:
                print(f"[FILE TYPE CHECK] Error fetching metadata for {srp_id}: {e}")

    # Safety net: any run_accession in gse_dict that we didn't classify
    # defaults to 'fastq' so existing behavior is preserved
    for gse_key in gse_dict.keys():
        for acc in gse_dict[gse_key]['run_accession']:
            if acc not in file_type_map:
                file_type_map[acc] = 'fastq'
                print(f"  {acc} -> FASTQ (default, no metadata match found)")

    # Summary
    bam_count = sum(1 for v in file_type_map.values() if v == 'bam')
    fastq_count = sum(1 for v in file_type_map.values() if v == 'fastq')
    print(f"\n[FILE TYPE CHECK] Summary: {fastq_count} FASTQ runs, {bam_count} BAM runs")

    return file_type_map


print("\n" + "=" * 80)
print("CLASSIFYING RUN ACCESSIONS BY ORIGINAL FILE TYPE")
print("=" * 80)
file_type_map = build_file_type_map(gse_dict)


#######################################################################
# BAM download and conversion helper functions
#######################################################################

def download_original_bam(accession, output_path):
    """
    Download the original BAM file for an SRA accession using prefetch --type all.
    Returns the path to the downloaded BAM file, or None if download failed.
    """
    print(f"[BAM DOWNLOAD] Attempting to download original BAM for {accession}")

    try:
        # prefetch with --type all downloads original submitted files (including BAMs)
        result = subprocess.Popen(
            ["prefetch", str(accession), "-O", str(output_path),
             "--max-size", "u", "--type", "all"],
            stdout=subprocess.PIPE, text=True
        )
        output, error = result.communicate()
        print(f'prefetch --type all outputs: {output}')
        print(f'prefetch --type all errors: {error}')

        # Search for the downloaded BAM file in the output directory
        bam_search_dir = os.path.join(str(output_path), str(accession))
        bam_files = glob.glob(os.path.join(bam_search_dir, '**', '*.bam'), recursive=True)

        if bam_files:
            print(f"[BAM DOWNLOAD] Found BAM file(s): {bam_files}")
            return bam_files[0]
        else:
            print(f"[BAM DOWNLOAD] No BAM files found in {bam_search_dir}")
            # Also check one level up in case directory structure differs
            bam_files = glob.glob(os.path.join(str(output_path), '**', '*.bam'), recursive=True)
            if bam_files:
                print(f"[BAM DOWNLOAD] Found BAM file(s) in broader search: {bam_files}")
                return bam_files[0]
            return None

    except Exception as e:
        print(f"[BAM DOWNLOAD] Failed to download BAM for {accession}: {e}")
        return None


def convert_bam_to_fastq(bam_path, accession, output_path, threads):
    """
    Convert a BAM file to FASTQ files using samtools.

    Produces:
      {accession}_1.fastq.gz (Read 1)
      {accession}_2.fastq.gz (Read 2)

    These match the naming convention that parallel-fastq-dump produces,
    so the existing renaming logic downstream works unchanged.
    """
    r1_path = os.path.join(str(output_path), f"{accession}_1.fastq.gz")
    r2_path = os.path.join(str(output_path), f"{accession}_2.fastq.gz")

    print(f"[BAM CONVERT] Converting {bam_path} to FASTQ files")
    print(f"[BAM CONVERT] Output R1: {r1_path}")
    print(f"[BAM CONVERT] Output R2: {r2_path}")

    try:
        # Step 1: Sort BAM by read name (required for proper paired-end extraction)
        sorted_bam = os.path.join(str(output_path), f"{accession}_namesorted.bam")
        print(f"[BAM CONVERT] Sorting BAM by read name...")

        sort_cmd = subprocess.Popen(
            ["samtools", "sort", "-n", "-@", str(threads),
             "-o", sorted_bam, str(bam_path)],
            stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True
        )
        output, error = sort_cmd.communicate()
        if sort_cmd.returncode != 0:
            print(f"[BAM CONVERT] samtools sort error: {error}")
            return False
        print(f"[BAM CONVERT] Sort complete.")

        # Step 2: Extract FASTQ from the name-sorted BAM
        print(f"[BAM CONVERT] Extracting FASTQ from sorted BAM...")

        r1_uncompressed = os.path.join(str(output_path), f"{accession}_1.fastq")
        r2_uncompressed = os.path.join(str(output_path), f"{accession}_2.fastq")

        fastq_cmd = subprocess.Popen(
            ["samtools", "fastq", "-@", str(threads),
             "-1", r1_uncompressed,
             "-2", r2_uncompressed,
             "-0", "/dev/null",
             "-s", "/dev/null",
             "-n",
             sorted_bam],
            stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True
        )
        output, error = fastq_cmd.communicate()
        if fastq_cmd.returncode != 0:
            print(f"[BAM CONVERT] samtools fastq error: {error}")
            return False
        print(f"[BAM CONVERT] FASTQ extraction complete.")

        # Step 3: Gzip the FASTQ files
        for fq_file in [r1_uncompressed, r2_uncompressed]:
            if os.path.exists(fq_file) and os.path.getsize(fq_file) > 0:
                print(f"[BAM CONVERT] Compressing {fq_file}...")
                gzip_cmd = subprocess.Popen(
                    ["gzip", "-f", fq_file],
                    stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True
                )
                gzip_cmd.communicate()
            else:
                print(f"[BAM CONVERT] Warning: {fq_file} is empty or missing")

        # Step 4: Clean up intermediate files
        if os.path.exists(sorted_bam):
            os.remove(sorted_bam)
            print(f"[BAM CONVERT] Cleaned up sorted BAM: {sorted_bam}")

        # Verify output files exist
        if os.path.exists(r1_path) and os.path.exists(r2_path):
            print(f"[BAM CONVERT] Successfully created FASTQ files from BAM")
            return True
        else:
            print(f"[BAM CONVERT] Warning: Expected output files not found after conversion")
            if os.path.exists(r1_path):
                print(f"[BAM CONVERT] Only R1 found - may be single-end data")
                return True
            return False

    except FileNotFoundError:
        print(f"[BAM CONVERT] ERROR: samtools not found. Please ensure samtools is installed.")
        print(f"[BAM CONVERT] Add samtools to the fastq_download_and_qc.yaml conda environment.")
        return False
    except Exception as e:
        print(f"[BAM CONVERT] ERROR during conversion: {e}")
        return False


#######################################################################
# Main download loop — existing FASTQ behavior preserved,
# BAM path only taken when metadata says the file IS a BAM.
#######################################################################

for key in gse_dict.keys():
    for accession in gse_dict[key]['run_accession']:
        print(f"\nProcessing sample {accession} from the BioProject {key}")

        sample_dir = os.path.join(str(output_dir), 'fastq', str(key), str(accession))

        # ------------------------------------------------------------------
        # BRANCH: Check file type map to decide download strategy
        # ------------------------------------------------------------------
        if file_type_map.get(accession, 'fastq') == 'bam':
            # ============================================================
            # BAM PATH: download original BAM, convert to FASTQ
            # ============================================================
            print(f"[BAM PATH] Metadata indicates {accession} was submitted as a BAM file")
            print(f"[BAM PATH] Downloading original BAM and converting to FASTQ...")

            # Ensure the sample output directory exists
            os.makedirs(sample_dir, exist_ok=True)

            bam_path = download_original_bam(accession, sample_dir)

            if bam_path:
                success = convert_bam_to_fastq(
                    bam_path, accession, sample_dir, computing_threads
                )
                if success:
                    print(f"[BAM PATH] Successfully converted BAM to FASTQ for {accession}")
                    # Clean up the original BAM to save disk space
                    try:
                        os.remove(bam_path)
                        print(f"[BAM PATH] Cleaned up original BAM: {bam_path}")
                    except Exception as e:
                        print(f"[BAM PATH] Could not remove BAM file: {e}")
                else:
                    print(f"[BAM PATH] WARNING: BAM to FASTQ conversion failed for {accession}")
                    print(f"[BAM PATH] This sample may need manual intervention.")
            else:
                print(f"[BAM PATH] WARNING: Could not download BAM for {accession}")
                print(f"[BAM PATH] This sample may need manual intervention.")

            # cd into sample dir to match existing behavior for the rename step
            os.chdir(sample_dir)

        else:
            # ============================================================
            # FASTQ PATH: existing behavior — prefetch + parallel-fastq-dump
            # ============================================================
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

        # ------------------------------------------------------------------
        # SHARED: Rename files to Illumina convention (runs for both paths)
        # ------------------------------------------------------------------
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
