#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os, sys

os.chdir(os.environ['HOME']+'/github/SRAscraper/tmp')

#%% Import the .acc accession firl downloaded from https://www.ncbi.nlm.nih.gov/labs/virus/vssi/#/virus?SeqType_s=Nucleotide&VirusLineage_ss=taxid:10239&HostLineage_ss=Homo%20sapiens%20(human),%20taxid:9606
import csv
ls_accessions = []
with open(os.environ['HOME']+'/github/SRAscraper/tmp/sequences_assembled.acc', 'r') as f:
          f_accessions = csv.reader(f, delimiter="\n")
          for row in f_accessions:
              ls_accessions.append(row[0])

len(ls_accessions)   
#%% Download the viral genomes using ncbi-datasets-cli

from multiprocessing.pool import ThreadPool
from math import ceil
from time import perf_counter
from os import cpu_count
import subprocess

NCPUS = cpu_count()

os.chdir(os.environ['HOME']+'/scPathoQuant/ref/tmp')

# Number of items in the list

def viral_scrape(accession):
    try:
        subprocess_6 = subprocess.Popen(
            ["datasets", "download", "virus", "genome", "accession", accession, '--filename', accession+'.zip'],
            stdout=subprocess.PIPE, text=True)
        output, error = subprocess_6.communicate()
        print(f'Outputs subprocess 6: {output}')
        print(f'Errors subprocess 6: {error}')
        subprocess_7 = subprocess.Popen(
            ['unzip', accession+'.zip', '-d', accession],
            stdout=subprocess.PIPE, text=True)
        output, error = subprocess_7.communicate()
        print(f'Errors subprocess 7: {error}')
        os.rename(os.environ['HOME']+'/scPathoQuant/ref/tmp/'+accession+'/ncbi_dataset/data/genomic.fna', os.environ['HOME']+'/scPathoQuant/ref/tmp/'+accession+'.fna')
        subprocess_8 = subprocess.Popen(
            ['rm -r '+accession+' '+accession+'.zip'],
            stdout=subprocess.PIPE, shell=True)
        output, error = subprocess_8.communicate()
        print(f'Errors subprocess 8: {error}')
    except:
        pass
            

#%%

with ThreadPool(NCPUS) as pool:
    chunksize = ceil(len(ls_accessions) / NCPUS)
    start = perf_counter()
    pool.map_async(viral_scrape, ls_accessions, chunksize=chunksize).get()
    end = perf_counter()
    print(f'Duration={end-start:.4f}s')


#%% out all the individual ref genomes toegther

subprocess_9 = subprocess.Popen(
    ['for file in *.fna; do cat \"$file\" >> viral_ref.fasta; done'],
    stdout=subprocess.PIPE, shell=True)
output, error = subprocess_9.communicate()
print(f'Errors: {error}')

#%%

subprocess_10 = subprocess.Popen(
    ['for file in *.fna; do rm file; done'],
stdout=subprocess.PIPE, shell=True)
output, error = subprocess_10.communicate()
print(f'Errors: {error}')


#%%

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
