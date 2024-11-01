#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os, sys

os.chdir(os.environ['HOME']+'/github/SRAscraper/tmp')

#%% Import in the dictionary with the metadata
import pickle

with open('dictionary_file.pkl', 'rb') as pkl_file:
     gse_dict = pickle.load(pkl_file)
     print('Dictionary loaded successfully')

     
#%% Run scPathoquant
import subprocess

from os import cpu_count

NCPUS = cpu_count()

cancer=("Bladder_cancer")
sample=("_S1_L001_")

subprocess_6 = subprocess.Popen(
    ["mkdir", "-p", os.environ['HOME']+"/scPathoQuant/results_"+cancer+"/pdf_bucket"],
    stdout=subprocess.PIPE, text=True)
output, error = subprocess_6.communicate()
print(f'Errors: {error}')

for key in gse_dict.keys():
    for accession in gse_dict[key]['SRR']:
        print(f"\nQuantifying viral reads in sample {accession} from the BioProject {key}")
        subprocess_7 = subprocess.Popen(
            ["mkdir", "-p", os.environ['HOME']+"/scPathoQuant/results_"+cancer+"/"+key+"/"+accession],
            stdout=subprocess.PIPE, text=True)
        output, error = subprocess_7.communicate()
        print(f'Errors: {error}')
        subprocess_8 = subprocess.Popen(
            ["scpathoquant", "-10x", os.environ['HOME']+"/SC/fastq/"+cancer+"/"+key+"/"+accession+"/"+accession+sample,
             "-op", os.environ['HOME']+"/scPathoQuant/results_"+cancer+"/"+key+"/"+accession,
             "-p", str(NCPUS), "-p2genome", os.environ['HOME']+"/scPathoQuant/ref"],
            stdout=subprocess.PIPE, text=True)
        output, error = subprocess_8.communicate()
        print(f'Errors: {error}')
        print(f'Outputs: {output}')
        print(f'Errors: {error}')
        subprocess_9 = subprocess.Popen(
            ["cd", os.environ['HOME']+"/scPathoQuant/results_"+cancer+"/"+key+"/"+accession],
            stdout=subprocess.PIPE, text=True)
        output, error = subprocess_9.communicate()
        print(f'Errors: {error}')
        subprocess_10 = subprocess.Popen(
            ['PDF_FILE=\"ls *.pdf\"; for f in $PDF_FILE; do cp \"$f\" '+os.environ['HOME']+'/scPathoQuant/results_'+cancer+'/pdf_bucket/'+accession+'_$f; done'],
            stdout=subprocess.PIPE, shell=True)
        output, error = subprocess_10.communicate()
        print(f'Errors: {error}')



#%% scPathoquant end
from Bio import SeqIO

fasta_sequences = SeqIO.to_dict(SeqIO.parse(os.environ['HOME']+'/scPathoQuant/ref/Ref-seq.fasta','fasta'))
fasta_sequences.keys()


#%% Something to gunzip the features file
import gzip
import csv


# read the tsv
for key in gse_dict.keys():
    for accession in gse_dict[key]['SRR']:
        NC_list = []
        with gzip.open(os.environ['HOME']+"/scPathoQuant/results_"+cancer+"/"+key+"/"+accession+"/filtered_feature_bc_matrix/features.tsv.gz", 'rt') as f:
            tsv_reader = csv.reader(f, delimiter="\t")
            for row in tsv_reader:
                NC_list.append(row[0])
            NC_list_final = [x for x in NC_list if x.startswith("NC_")]

    
            # add rows where there are not rows
            for virus in fasta_sequences.keys():
                if virus in NC_list_final:
                    print(f'{virus} reads found in {accession}')
                else:
                    with gzip.open(os.environ['HOME']+"/scPathoQuant/results_"+cancer+"/"+key+"/"+accession+"/filtered_feature_bc_matrix/features.tsv.gz", 'at') as f:
                        tsv_writer = csv.writer(f, delimiter="\t", lineterminator="\n")
                        tsv_writer.writerow([virus, virus, 'Gene Expression'])\
                
                
        with gzip.open(os.environ['HOME']+"/scPathoQuant/results_"+cancer+"/"+key+"/"+accession+"/filtered_feature_bc_matrix/features.tsv.gz", 'rt') as f:
            feature_tsv = []
            tsv_reader = csv.reader(f, delimiter="\t")
            for row in tsv_reader:
                feature_tsv.append(row[0])
                
                
            feature_tsv_len = len(feature_tsv)
        with gzip.open(os.environ['HOME']+"/scPathoQuant/results_"+cancer+"/"+key+"/"+accession+"/filtered_feature_bc_matrix/matrix.mtx.gz") as file:
            lines = file.readlines()
            tmp = str(feature_tsv_len) + ' ' + str(lines[3 - 1]).split(' ', 2)[1] + ' ' + str(len(lines) - 3) + '\n'
            lines[3 - 1] = bytes(tmp, 'utf-8')
            
            
        with gzip.open(os.environ['HOME']+"/scPathoQuant/results_"+cancer+"/"+key+"/"+accession+"/filtered_feature_bc_matrix/matrix.mtx.gz", "wb") as file:
            for line in lines:
                file.write(line)


#%% End file

sys.exit()
