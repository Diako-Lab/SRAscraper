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

     
# Go get some files
import subprocess

for key in gse_dict.keys():
    for accession in gse_dict[key]['run_accession']:
        print(f"\nProcessing sample {accession} from the BioProject {key}")
        subprocess_1 = subprocess.Popen(
            ["parallel-fastq-dump", "--sra-id", accession, "--threads", str(computing_threads), "--outdir", 
             output_dir+'/fastq/'+key+'/'+accession, "--split-spot", "--split-files", "--gzip"], stdout=subprocess.PIPE, text=True)
        output, error = subprocess_1.communicate()
        print(f'Outputs: {output}')
        print(f'Errors: {error}')
        print(f"\nRenamming {accession} fastqs")
        try:
             os.rename(output_dir + '/fastq/' + key + '/' + accession + '/' + accession + '_3.fastq.gz',
                    output_dir + '/fastq/' + key + '/' + accession + '/' + accession + '_S1_L001_I1_001.fastq.gz')
        except FileNotFoundError:
            print(f"\nIndex {accession}_3.fastq.gz file not found in GEO repo. Continuing with normal R1 R2 renaming.")
        except OSError as e:
            print("Error renaming file:", e)
        try:
            os.rename(output_dir + '/fastq/' + key + '/' + accession + '/' + accession + '_1.fastq.gz',
                    output_dir + '/fastq/' + key + '/' + accession + '/' + accession + '_S1_L001_R1_001.fastq.gz')
            os.rename(output_dir + '/fastq/' + key + '/' + accession + '/' + accession + '_2.fastq.gz',
                    output_dir + '/fastq/' + key + '/' + accession + '/' + accession + '_S1_L001_R2_001.fastq.gz')
        except FileNotFoundError:
            print(f"\nIndex {accession}_1/2.fastq.gz files not found in GEO repo.")
        except OSError as e:
            print("Error renaming file:", e)



# QC Section
# Go check some files with fastqc

for key in gse_dict.keys():
    for accession in gse_dict[key]['run_accession']:
        subprocess_2 = subprocess.Popen(
            ["fastqc", "-t", str(computing_threads), "-o", os.path.join(output_dir, 'QC'), os.path.join(output_dir, 'fastq', key, accession, accession+'_S1_L001_R1_001.fastq.gz'), os.path.join(output_dir, 'fastq',  key, accession, accession+'_S2_L001_R1_001.fastq.gz')], stdout=subprocess.PIPE, text=True)
        output, error = subprocess_2.communicate()
        print(f'Outputs: {output}')
        print(f'Errors: {error}')

# Pull all the repoorts together with multiqc

subprocess_3 = subprocess.Popen(
    ["multiqc", "-o", os.path.join(output_dir, 'QC'), os.path.join(output_dir, 'QC')], stdout=subprocess.PIPE, text=True)
output, error = subprocess_3.communicate()
print(f'Outputs: {output}')
print(f'Errors: {error}')


# End file

sys.exit()
