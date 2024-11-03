#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os, sys

output_dir = snakemake.params.output_dir
NCBI_search_txt = snakemake.params.NCBI_search_txt
computing_threads = snakemake.params.computing_threads

os.chdir(output_dir)

#%% Srape out the ftp links from the NCBI search results we will use this later
import pandas as pd
pd.options.display.max_colwidth = 10000
# Read in the NCBI search result file
X = pd.read_fwf(NCBI_search_txt, header=None)


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


#%% Multiprocessed code

from multiprocessing.pool import ThreadPool
from math import ceil
from time import perf_counter

ftp_list_input = ftp_list

# Number of items in the list
N = len(ftp_list_input)

with ThreadPool(computing_threads) as pool:
    chunksize = ceil(len(ftp_list_input) / computing_threads)
    start = perf_counter()
    results = list(pool.map_async(url_ok, ftp_list_input, chunksize=chunksize).get())
    end = perf_counter()
    assert len(results) == N
    print(f'Duration={end-start:.4f}s')


def countX(lst, x):
    count = 0
    for ele in lst:
        if (ele == x):
            count = count + 1
    return count

print(f'''Of the initial imported {len(ftp_links)} datasets provided from the NCBI search,
      {len(ftp_list)} point to GEO FTP pages.
      Of these {len(ftp_list)} FTP pages {countX(results, True)} returned 200 GET statuses for having GEO *.soft files.''')


from itertools import compress
ftp_list = list(compress(ftp_list, results))
len(ftp_list)

# Go fetch all the valid GEO ftp links and get the .soft files
# These will be needed to fetch all the GEO ftp files 
import wget
import GEOparse

# These are needed to get the GSE files to SRR file names
import subprocess
subprocess.run(["pysradb", "--version"])

from pysradb.search import GeoSearch
# Make a big dictionary that contains all the .soft file information
gse_dict = {}
for url in ftp_list:
    try:
        wget_this = url + 'soft/' + url[:-1].split('/', -1)[-1] + '_family.soft.gz'
        filename = wget.download(wget_this)
        gse = GEOparse.get_GEO(filepath=filename)
        gse_dict[url[:-1].split('/', -1)[-1]] = gse.phenotype_data
        instance = GeoSearch(2, 200, geo_query=url[:-1].split('/', -1)[-1])
        instance.search()
        df = instance.get_df()
        data_pysrad = {'SRR': df.run_1_accession.tolist(), 
                'geo_accession': df.sample_alias.tolist(),
                'Rep': list(range(1, len(df.run_1_accession)+1))}
        df_pysrad = pd.DataFrame(data_pysrad)
        gse_dict[url[:-1].split('/', -1)[-1]] = gse_dict[url[:-1].split('/', -1)[-1]].set_index('geo_accession').merge(df_pysrad, right_on='geo_accession', left_index=True)
    except:
        pass

# Call out the number of bio samples that are available from the .soft files before filtering
num_rows = []
unique_columns = []
for key in gse_dict.keys():
    num_rows.append(len(gse_dict[key]))
    unique_columns = [*unique_columns, *gse_dict[key].columns]

print(f'''From the {countX(results, True)} FTP pages that returned 200 GET 
      statuses for having GEO *.soft files. 
      There are an available {sum(num_rows)} samples that are available for download.''')

# Save the file This would be the stopping point for the first snakemake rule
import pickle

with open('dictionary_file.pkl', 'wb') as pkl_file:
    pickle.dump(gse_dict, pkl_file)
    print('Dictionary saved successfully.')
    

# End file
import sys

sys.exit()
