#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os, sys

output_dir = snakemake.params.output_dir
computing_threads = snakemake.params.computing_threads

os.chdir(output_dir)


#%% Import in the dictionary with the metadata
import pickle

with open('dictionary_file.pkl', 'rb') as pkl_file:
     gse_dict = pickle.load(pkl_file)
     print('Dictionary loaded successfully')

     
#%% Go get some files
import subprocess

for key in gse_dict.keys():
    for accession in gse_dict[key]['SRR']:
        print(f"\nProcessing sample {accession} from the BioProject {key}")
        subprocess_1 = subprocess.Popen(

            ["parallel-fastq-dump", "--sra-id", accession, "--threads", computing_threads, "--outdir", 
             output_dir + '/fastq/'+key+'/'+accession, "--split-spot", "--split-files", "--gzip"], stdout=subprocess.PIPE, text=True)

        output, error = subprocess_1.communicate()
        print(f'Outputs: {output}')
        print(f'Errors: {error}')
        print(f"\nRenamming {accession} fastqs")

        try:
            os.rename(output_dir + '/fastq/' + key + '/' + accession + '/' + accession + '_1.fastq.gz',
                    output_dir + '/fastq/' + key + '/' + accession + '/' + accession + '_S1_L001_R1_001.fastq.gz')
            os.rename(output_dir + '/fastq/' + key + '/' + accession + '/' + accession + '_2.fastq.gz',
                    output_dir + '/fastq/' + key + '/' + accession + '/' + accession + '_S2_L001_R1_001.fastq.gz')
        except FileNotFoundError:
            print("File not found.")
        except OSError as e:
            print("Error renaming file:", e)



#%% End file

sys.exit()
