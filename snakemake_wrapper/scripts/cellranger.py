#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os, sys

os.chdir(os.environ['HOME']+'/github/SRAscraper/tmp')

#%% Import in the dictionary with the metadata
import pickle

with open('dictionary_file.pkl', 'rb') as pkl_file:
     gse_dict = pickle.load(pkl_file)
     print('Dictionary loaded successfully')


#%% Yehaw, cellranger time
import subprocess

from os import cpu_count

NCPUS = cpu_count()

cancer=("Bladder_cancer")
sample=("_S1_L001_")

for key in gse_dict.keys():
    for accession in gse_dict[key]['SRR']:
        os.chdir(os.environ['HOME']+'/SC/fastq/'+cancer+'/'+key+'/'+accession+'/')
        subprocess_5 = subprocess.Popen(
            ["cellranger", "count", "--id="+accession+sample, 
             "--fastqs="+os.environ['HOME']+'/SC/fastq/'+cancer+'/'+key+'/'+accession,
             "--sample="+accession,
             "--transcriptome="+os.environ['HOME']+'/SC/ref/GRCh38',
             "--localcores="+str(NCPUS),
             "--create-bam=true",
             "--chemistry=auto"], stdout=subprocess.PIPE, text=True)
        output, error = subprocess_5.communicate()
        print(f'Outputs: {output}')
        print(f'Errors: {error}')


#%% End file

sys.exit()
