#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 17 09:35:27 2024

@author: jlehle
"""
# %%


import pandas as pd
from bs4 import BeautifulSoup
import xml.etree.ElementTree as ET
from random import randint
from random import seed
import time
import json
import csv
import re




import tarfile

#%%
import os, sys
#os.chdir(os.environ['HOME']+'/SC/fastq')
os.chdir(os.environ['HOME']+'/OmniScrape/tmp')
#%%
#######
# THIS IS WHERE THINGs STARTED WORKING #
#######

#%% Srape out the ftp links from the SRA or GEO search results we will use this later

import pandas as pd
pd.options.display.max_colwidth = 1000
#X = pd.read_fwf('/master/jlehle/OmniScrape/tmp/gds_result.txt', header=None)
#X = pd.read_fwf('/master/jlehle/OmniScrape/tmp/Smart-seq-GEO.txt', header=None)
X = pd.read_fwf('/master/jlehle/OmniScrape/tmp/kidney_tumor_SRA.txt', header=None)
X[0]

ftp_links = X[X[0].str.startswith('FTP download')]
ftp_links = ftp_links.reset_index(drop=True)
ftp_links = ftp_links.loc[0:,0]
ftp_links


# %% scrape the ftp from the gds file
ftp_list = []
# loop through the rows using iterrows()
for index, row in ftp_links.items():
    try:
        tmp = row.split(') ', -1)[1]
        ftp_list.append(tmp)
    except:
        pass


# %% check the ftp links
# import module
import requests
# create a function
# pass the url

def url_ok(foo_url):
    foo_url = 'https' + foo_url[3:] + 'soft/'
    # pass the url into
    # request.hear
    response = requests.head(foo_url)
    # check the status code
    if response.status_code == 200:
        return True
    else:
        return False


#%%
# Multiprocessed code

from multiprocessing.pool import ThreadPool
from math import ceil
from time import perf_counter
from os import cpu_count


ftp_list_input = ftp_list

# Number of items in the list
N = len(ftp_list_input)

try: 
    maunal_NCPUS
    print(f'Creating pool of {manual_NCPUS} threads')
except NameError:
    # If no manual number set use half of available CPUs (minimum of 2)
    if (NCPUS := cpu_count()):
        NCPUS = max(NCPUS//2, 2)
        print(f'Creating pool of {NCPUS} threads')
    else:
        NCPUS = 2
        print(f'Creating pool of {NCPUS} threads')

with ThreadPool(NCPUS) as pool:
    chunksize = ceil(len(ftp_list_input) / NCPUS)
    start = perf_counter()
    results = list(pool.map_async(url_ok, ftp_list_input, chunksize=chunksize).get())
    end = perf_counter()
    assert len(results) == N
    print(f'Duration={end-start:.4f}s')


#%%
def countX(lst, x):
    count = 0
    for ele in lst:
        if (ele == x):
            count = count + 1
    return count


#%%
print(f'''Of the initiaimport wgetl {len(ftp_links)} datasets provided from the GEO search,
      {len(ftp_list)} point to GEO FTP pages.
      Of these {len(ftp_list)} FTP pages {countX(results, True)} returned 200 GET statuses for having GEO *.soft files.''')

#%%
from itertools import compress
ftp_list = list(compress(ftp_list, results))
len(ftp_list)
#%% Go fetch all the valid GEO ftp links and get the .soft files
import wget
import GEOparse
# Make a big dictionary that contains all the .soft file information
gse_dict = {}
for url in ftp_list:
    wget_this = url + 'soft/' + url[:-1].split('/', -1)[-1] + '_family.soft.gz'
    filename = wget.download(wget_this)
    gse = GEOparse.get_GEO(filepath=filename)
    gse_dict[url[:-1].split('/', -1)[-1]] = gse.phenotype_data


#%% Call out the number of bio samples that are available from the .soft files before filtering
num_rows = []
unique_columns = []
for key in gse_dict.keys():
    num_rows.append(len(gse_dict[key]))
    unique_columns = [*unique_columns, *gse_dict[key].columns]

#%%
print(f'''From the {countX(results, True)} FTP pages that returned 200 GET 
      statuses for having GEO *.soft files. 
      There are an available {sum(num_rows)} samples that are available for download.''')

######

###STOP HERE??

#%% WORKING WORKING

gse_dict.keys()
[gse_dict['GSE73121']['title']]

test = pd.DataFrame.from_dict(gse_dict['GSE73121'])
test
unique_columns

df_chunk = []
df_final = []
for key in gse_dit.keys():
        for index, item in enumerate(np.unique(unique_columns)):
            try:
                data = [gse_dict[key][item]]
                df = pd.DataFrame(data, columns=[item])
                df_chunk = pd.concat([df_chunk.reset_index(drop=True), df.reset_index(drop=True)], axis=1)
            else:
                if item not in gse_dict[key].columns:
                    ####df[item] = np.nan



data = {'Project': np.repeat('Espophageal Cancer', [128]).tolist(),
        'Sample': gse.phenotype_data.title.tolist(),
        'GSM': gse.phenotype_data.geo_accession.tolist(),
        'Condition': gse.phenotype_data.iloc[0:,11].tolist(),
        'Source': gse.phenotype_data.source_name_ch1.tolist(),
        'CD45': gse.phenotype_data.iloc[0:,12].tolist(),
        }
df_geo = pd.DataFrame(data)



#Start filtering the 25,202 samples

#Figure out some way to write everything you have collected so far into table this will be the first target file of the script. 
I will need to have GSE the 



GSE_dic_final = dict()
for key in gse_dict.keys():
    for index, item in enumerate(np.unique(unique_columns)):
        if item not in df.columns:
            df[item] = np.nan
    GSE_dic = {row: df}
    GSE_dic_final = GSE_dic_final | GSE_dic






#%% This section is correct as well It is where I pull the SRA informtion from 
#It has to start with the SRA or GEO search results

GSE = X[X[0].str.startswith('Series		Accession: GSE')]
GSE = GSE.reset_index(drop=True)
GSE = GSE.iloc[0:,0]
GSE

#%% scrape the GSE from the gds file
GSE_list = []
# loop through the rows using iterrows()
for index, row in GSE.items():
    try:
        tmp = row.split('Series\t\tAccession: ', 1)[1]
        tmp = tmp.split('\t', 1)[0]
        GSE_list.append(tmp)
    except:
        pass

#%%
import subprocess
subprocess.run(["pysradb", "--version"])
from pysradb.search import GeoSearch
import numpy as np

unique_columns = []

for index, row in enumerate(GSE_list):
    print(row)
    instance = GeoSearch(2, 1000, geo_query=row)
    instance.search()
    df = instance.get_df()
    unique_columns = [*unique_columns, *df.columns.to_list()]
    print(unique_columns)
    

#%%  
GSE_dic_final = dict()
for index, row in enumerate(GSE_list):  
    instance = GeoSearch(2, 1000, geo_query=row)
    instance.search()
    df = instance.get_df()
    for index, item in enumerate(np.unique(unique_columns)):
        if item not in df.columns:
            df[item] = np.nan
    GSE_dic = {row: df}
    GSE_dic_final = GSE_dic_final | GSE_dic
    
    

#%%

sample = []
organism = []
library_se = []
library_st = []
gse_dict_copy = gse_dict
for key in gse_dict_copy.keys():
    gse = GEOparse.get_GEO(filepath= key+"_family.soft.gz")
    for index, gsm in gse.phenotype_data['geo_accession'].items():
        sample.append(gse.phenotype_data.loc[gsm]['geo_accession'])
        organism.append(gse.phenotype_data.loc[gsm]['organism_ch1'])
        library_se.append(gse.phenotype_data.loc[gsm]['library_selection'])
        library_st.append(gse.phenotype_data.loc[gsm]['library_strategy'])
    
#%%  
# function to get unique values
def unique(list1):
    unique_list = pd.Series(list1).drop_duplicates().tolist()
    for x in unique_list:
        print(x)
#%%
organism
unique(organism)
unique(library_se)
unique(library_st)
type(organims)






#%%
##### Leftovers from setting up the GSE_dic section
XXXXXXXXXXXXX
XXXXXXXXXXXX





GSE_dic_final['GSE182632'].iloc[0]

len(unique_columns)


    ######## Stoped here!!!! put in starting dataframe
print(df.iloc[0])
#%% Esophageal cancer datasets
import subprocess
subprocess.run(["pysradb", "--version"])

from pysradb.search import GeoSearch
instance = GeoSearch(2, 200, geo_query='GSE173343')
instance.search()
df = instance.get_df()
type(df)
GSE_dic = {'GSE173343': df}
GSE_dic['GSE173343'].run_1_accession
tmp = df.columns.to_list()
tmp[0]
df.columns.to_list()
print(df.experiment_title)

data_pysrad = {'SRR': df.run_1_accession.tolist(), 
        'GSM': df.sample_alias.tolist(),
        'Rep': list(range(1, len(df.run_1_accession)+1))}
df_pysrad = pd.DataFrame(data_pysrad)
df_pysrad
#%%
import wget
import GEOparse

wget_this = 'ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE160nnn/GSE160269/soft/GSE160269_family.soft.gz'
filename = wget.download(wget_this)
gse = GEOparse.get_GEO(filepath=filename)
#%%
adata_pp
#%%
data = {'Project': np.repeat('Espophageal Cancer', [128]).tolist(),
        'Sample': gse.phenotype_data.title.tolist(),
        'GSM': gse.phenotype_data.geo_accession.tolist(),
        'Condition': gse.phenotype_data.iloc[0:,11].tolist(),
        'Source': gse.phenotype_data.source_name_ch1.tolist(),
        'CD45': gse.phenotype_data.iloc[0:,12].tolist(),
        }
df_geo = pd.DataFrame(data)

df = df_pysrad.set_index('GSM').join(df_geo.set_index('GSM'))
df

#Path to the files. Look how the file ends with the id beginning of the file and the only thing missing is the sample id value that is unique for each sample.
#file_base = '/master/jlehle/SC/fastq/Esophagus_cancer/'
file_base = '/master/jlehle/scPathoQuant/results_Esophagus_cancer/'

df








#%%
#######
# THIS IS WHERE I PUT THE CODE THAT I'M WORKING ON #
#######

#%%
#Mohadeseh needs more samples so I'm gonna set that up really fast

wget_this = 'ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE160nnn/GSE160269/soft/GSE160269_family.soft.gz'
filename = wget.download(wget_this)
gse = GEOparse.get_GEO(filepath=filename)
gse.phenotype_data
SRA_file = pd.read_csv("test.tsv", sep="\t")
SRA_file.loc[0]
SRA_file.run_accession
print(SRA_file.experiment_title)
SRA_file.biosample

data = {'SRR': SRA_file.run_accession,
        'GSM': SRA_file.experiment_title}


#%% Bladder cancer datasets
import wget
import GEOparse

wget_this = 'ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE222nnn/GSE222315/soft/GSE222315_family.soft.gz'
filename = wget.download(wget_this)
gse = GEOparse.get_GEO(filepath=filename)

data = {'Sample': ["P1_BCa_scRNAseq", "P1_BCa_scRNAseq", "P1_BCa_scRNAseq", "P1_BCa_scRNAseq", \
"P2_NAT_scRNAseq", "P2_NAT_scRNAseq", "P2_NAT_scRNAseq", "P2_NAT_scRNAseq", \
"P2_BCa_scRNAseq", "P2_BCa_scRNAseq", "P2_BCa_scRNAseq", "P2_BCa_scRNAseq", \
"P3_NAT_scRNAseq", "P3_NAT_scRNAseq", "P3_NAT_scRNAseq", "P3_NAT_scRNAseq", \
"P3_BCa_scRNAseq", "P3_BCa_scRNAseq", "P3_BCa_scRNAseq", "P3_BCa_scRNAseq", \
"P4_NAT_scRNAseq", "P4_NAT_scRNAseq", "P4_NAT_scRNAseq", "P4_NAT_scRNAseq", \
"P4_BCa_scRNAseq", "P4_BCa_scRNAseq", "P4_BCa_scRNAseq", "P4_BCa_scRNAseq", \
"P5_NAT_scRNAseq", "P5_NAT_scRNAseq", "P5_NAT_scRNAseq", "P5_NAT_scRNAseq", \
"P5_BCa_scRNAseq", "P5_BCa_scRNAseq", "P5_BCa_scRNAseq", "P5_BCa_scRNAseq", \
"P6_BCa_scRNAseq", \
"P7_BCa_scRNAseq", \
"P8_BCa_scRNAseq", \
"P9_BCa_scRNAseq"],
        'Condition': np.repeat(gse.phenotype_data.source_name_ch1, [4, 4, 4, 4, 4, 4, 4, 4, 4, 1, 1, 1, 1]).tolist(),
        'Source': np.repeat('Bladder', [40]).tolist(),
        'Cancer Stage': np.repeat(gse.phenotype_data['characteristics_ch1.3.tumor stage'],  [4, 4, 4, 4, 4, 4, 4, 4, 4, 1, 1, 1, 1]).tolist(),
        'Sex': np.repeat(gse.phenotype_data['characteristics_ch1.1.Sex'], [4, 4, 4, 4, 4, 4, 4, 4, 4, 1, 1, 1, 1]).tolist(),
        'Age': np.repeat(gse.phenotype_data['characteristics_ch1.2.age'], [4, 4, 4, 4, 4, 4, 4, 4, 4, 1, 1, 1, 1]).tolist()       
        }
df = pd.DataFrame(data)
#%%
len(data['Sample'])
len(data['Source'])
df = pd.DataFrame(data)
df.Sample.tolist()
#%%    
    
    
SRA_file = pd.read_csv("test.tsv", sep="\t")
SRA_file.loc[0]
SRA_file.run_accession
print(SRA_file.experiment_title)
SRA_file.biosample

data = {'SRR': SRA_file.run_accession,
        'GSM': SRA_file.experiment_title}
#%%
new_string_final = []
for index, string in data['GSM'].items():
    new_string = string.split(' 10X', 1)[0]
    new_string_final.append(new_string)
#%%    
data['Sample'] = new_string_final    
data['SRR'].to_string()
    
gse_dict[url[:-1].split('/', -1)[-1]] = gse.phenotype_data

#%%
#the sample
len(gse.phenotype_data['geo_accession'])
gse.phenotype_data.loc['GSM3327702']['organism_ch1']
#organism is
gse.phenotype_data['organism_ch1']
library info
gse.phenotype_data['library_selection']
gse.phenotype_data['library_strategy']


# %% *****
gse_dict.keys()
gse = GEOparse.get_GEO(filepath="GSE118389_family.soft.gz")
type(gse.phenotype_data)
gse.phenotype_data
gse.gpls
gse.gsms
gse.gsms['GSM4909314'].show_metadata()
gse.gsms['GSM6129415'].get_metadata_attribute('Sample_title')
gse.gsms['GSM4909314'].download_SRA(
    'jlehle@txbiomed.org', directory=GSE_UID, nproc=1,)
gse.gsms['GSM4909314'].download_supplementary_files()

pd.set_option('display.max_columns', None)
print(gse.phenotype_data)
print(gse.phenotype_data['geo_accession'])

for accession in gse.phenotype_data['geo_accession']; :
    gse.gsms[accession].download_SRA(
        'jlehle@txbiomed.org', directory=gse.phenotype_data['series_id'], nproc=1)
print(gse.phenotype_data['series_id'])


# %%
#######
# THIS IS WHERE I PUT THE CODE THAT FAILED ME#
#######

# %%
def extract(tar_url, extract_path='.'):
    print(tar_url)
    tar = tarfile.open(tar_url, 'r')
    for item in tar:
        tar.extract(item, extract_path)
        if item.name.find(".tgz") != -1 or item.name.find(".tar") != -1:
            extract(item.name, "./" + item.name[:item.name.rfind('/')])


try:

    extract(sys.argv[1] + '.tgz')
    print('Done.')
except:
    name = os.path.basename(sys.argv[0])
    print(name[:name.rfind('.')], '<filename>')
# %%
extract('GSE252723_family.xml.tgz')

# %%
# %%
# Reading the data inside the xml
# file to a variable under the name
# data
# %%
with open('GSE252723_family.xml', 'r') as f:
    data = f.read()

# %%
# Passing the stored data inside
# the beautifulsoup parser, storing
# the returned object
Bs_data = BeautifulSoup(data, "xml")
print(Bs_data)
# %%
xmlTree = ET.parse('GSE252723_family.xml')
print(xmlTree)
tags = {elem.tag for elem in xmlTree.iter()}
tags
# %%
# Finding all instances of tag
# `unique`
b_unique = Bs_data.find_all('Data-Table')
b_unique
# %%
headers = {
    'authority': 'screen-beta-api.wenglab.org',
    'accept': 'application/json',
    'accept-language': 'en-US,en;q=0.9',
    'content-type': 'application/json',
    'origin': 'https://screen.wenglab.org',
    'referer': 'https://screen.wenglab.org/',
    'sec-ch-ua': '\".Not/A)Brand\";v=\"99\", \"Google Chrome\";v=\"103\", \"Chromium\";v=\"103\"',
    'sec-ch-ua-mobile': '?1',
    'sec-ch-ua-platform': '\"Android\"',
    'sec-fetch-dest': 'empty',
    'sec-fetch-mode': 'cors',
    'sec-fetch-site': 'same-site',
    'user-agent': 'Mozilla/5.0 (Linux; Android 6.0; Nexus 5 Build/MRA58N) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/103.0.0.0 Mobile Safari/537.36'
}
# %%
url = 'https://www.ncbi.nlm.nih.gov/gds/?term=200264681'
# %%
url_output = requests.get(url)
url_output.headers['Content-Type']
url_output.text
geo_json = url_output.json()

pd.readhttps: // www.ncbi.nlm.nih.gov/gds /?term = 200264681
type(X)
X
X.loc[0, 0].startswith('1.')


fname = 'guppy-0.1.10.tar.gz'
url = 'https://pypi.python.org/packages/source/g/guppy/' + fname
r = requests.get(url)
open(fname, 'wb').write(r.content)
