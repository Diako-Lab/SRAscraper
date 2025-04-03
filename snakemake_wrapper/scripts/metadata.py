#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os, sys

metadata_dir = os.path.join(snakemake.params.output_dir, 'metadata')
NCBI_search_txt = snakemake.params.NCBI_search_txt
computing_threads = snakemake.params.computing_threads

isExist = os.path.exists(metadata_dir)

if not isExist:
    os.mkdir(metadata_dir)


os.chdir(metadata_dir)

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
import time
# create a function
# pass the url

def url_ok(foo_url):
    transformed_url = 'https' + foo_url[3:] + 'soft/'
    # pass the url into
    # request.hear
    max_retries=10 
    retry_delay=1
    for attempt in range(max_retries):
        try:
            response = requests.head(transformed_url, timeout=10)
            if response.status_code == 200:
                return True
            print(f"Attempt {attempt+1}: Got status {response.status_code} for {transformed_url}")
        except Exception as e:
            print(f"Attempt {attempt+1} failed for {transformed_url}: {str(e)}")
        # If not 200 and not last attempt, wait and retry
        if attempt < max_retries - 1:
            time.sleep(retry_delay)
    
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

def download_with_retry(url, max_retries=10, retry_delay=1):
    """Download a file with retries on failure"""
    last_exception = None
    for attempt in range(max_retries):
        try:
            filename = wget.download(url)
            return filename
        except (HTTPError, URLError, Exception) as e:
            last_exception = e
            print(f"Attempt {attempt + 1} failed for {url}: {str(e)}")
            if attempt < max_retries - 1:
                time.sleep(retry_delay)
    
    print(f"All {max_retries} attempts failed for {url}")
    raise last_exception if last_exception else Exception("Download failed")
# Go fetch all the valid GEO ftp links and get the .soft files
# These will be needed to fetch all the GEO ftp files 
import wget
import GEOparse

# These are needed to get the GSE files to SRR file names
import subprocess
subprocess.run(["pysradb", "--version"])

from pysradb.search import GeoSearch
from pysradb.sraweb import SRAweb
sradb = SRAweb()
# Make a big dictionary that contains all the .soft file information
gse_dict = {}
for url in ftp_list:
    print(url)
    try:
        wget_this = url + 'soft/' + url[:-1].split('/', -1)[-1] + '_family.soft.gz'
        filename = download_with_retry(wget_this)
        gse = GEOparse.get_GEO(filepath=filename)
        gse_df = gse.phenotype_data
        try:
            columns_to_select = ['title', 'geo_accession', 'extract_protocol_ch1', 'description, data_processing', 'platform_id',	'contact_name', 'contact_department', 'contact_institute', 'contact_address', 'contact_city', 'contact_state', 'contact_zip/postal_code', 'contact_country']
            gse_df_filtered =  gse_df.drop(columns=gse_df.columns.difference(columns_to_select))
            instance = GeoSearch(2, 200, geo_query=url[:-1].split('/', -1)[-1])
            instance.search()
            instance_df = instance.get_df()
            instance_df
            srp_dfs = []
            for srp_id in instance_df.study_accession.unique():
                srp_df = sradb.sra_metadata(srp_id, detailed = True)
                srp_dfs.append(srp_df)
               
            if len(srp_dfs) > 1:
                final_srp_df = pd.concat(srp_dfs, ignore_index=True)
                #This might get you in trouble if the columns are different between srp look back here to fix that
            else:
                final_srp_df = srp_dfs[0]
            #This removes everything except the sample alias and run_1_accession
            instance_df.drop(columns=instance_df.columns.difference(['run_1_accession', 'sample_alias']), inplace=True)
            #Then rename the columns
            instance_df = instance_df.rename(columns={'run_1_accession': 'run_accession', 'sample_alias' : 'geo_accession'})
            tmp = pd.merge(final_srp_df, instance_df, on='run_accession', how='outer')
            final_srp_df = pd.merge(final_srp_df, instance_df, on='run_accession', how='outer')
            gse_dict[url[:-1].split('/', -1)[-1]] = gse_df_filtered.merge(final_srp_df, right_on='geo_accession', left_index=True, how='outer')
            gse_dict[url[:-1].split('/', -1)[-1]].drop(columns=['geo_accession_x', 'geo_accession_y'], inplace=True)
        except:
            gse_df = gse.phenotype_data
            gse_dict[url[:-1].split('/', -1)[-1]] = gse_df
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
