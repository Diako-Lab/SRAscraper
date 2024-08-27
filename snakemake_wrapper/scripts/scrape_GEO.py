#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 17 09:35:27 2024

@author: jlehle
"""
import os, sys

os.chdir(os.environ['HOME']+'/OmniScrape/tmp')

#%% Srape out the ftp links from the SRA or GEO search results we will use this later
import pandas as pd
pd.options.display.max_colwidth = 10000
# Reasign X to config once this is working
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
#%%
import yaml

with open('dictionary_file.yml', 'w') as yaml_file:
     yaml.dump(gse_dict, stream=yaml_file, default_flow_style=False)