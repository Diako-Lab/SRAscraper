#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os, sys

os.chdir(os.environ['HOME']+'/github/SRAscraper/tmp')

#%% Import in the dictionary with the metadata
import pickle

with open('dictionary_file.pkl', 'rb') as pkl_file:
     gse_dict = pickle.load(pkl_file)
     print('Dictionary loaded successfully')

     
#%% Go get some files
import subprocess

from os import cpu_count

NCPUS = cpu_count()

cancer=("Bladder_cancer")

for key in gse_dict.keys():
    for accession in gse_dict[key]['SRR']:
        print(f"\nProcessing sample {accession} from the BioProject {key}")
        subprocess_1 = subprocess.Popen(
            ["parallel-fastq-dump", "--sra-id", accession, "--threads", str(NCPUS), "--outdir", 
             os.environ['HOME']+'/SC/fastq/'+cancer+'/'+key+'/'+accession,
             "--split-spot", "--split-files", "--gzip"], stdout=subprocess.PIPE, text=True)
        output, error = subprocess_1.communicate()
        print(f'Outputs: {output}')
        print(f'Errors: {error}')
        print(f"\nRenamming {accession} fastqs")
        subprocess_2 = subprocess.Popen(
            ["mv", os.environ['HOME']+'/SC/fastq/'+cancer+'/'+key+'/'+accession+'/'+accession+'_1.fastq.gz',
             os.environ['HOME']+'/SC/fastq/'+cancer+'/'+key+'/'+accession+'/'+accession+'_S1_L001_R1_001.fastq.gz'], 
            stdout=subprocess.PIPE, text=True)
        output, error = subprocess_2.communicate()
        print(f'Errors: {error}')
        subprocess_3 = subprocess.Popen(
            ["mv", os.environ['HOME']+'/SC/fastq/'+cancer+'/'+key+'/'+accession+'/'+accession+'_2.fastq.gz',
              os.environ['HOME']+'/SC/fastq/'+cancer+'/'+key+'/'+accession+'/'+accession+'_S1_L001_R2_001.fastq.gz'], 
            stdout=subprocess.PIPE, text=True)
        output, error = subprocess_3.communicate()
        print(f'Errors: {error}')   

#%% End file

sys.exit()
