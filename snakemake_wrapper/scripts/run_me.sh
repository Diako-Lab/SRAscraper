#!/bin/bash

#PBS -N Monster_Food
#PBS -l nodes=1:ppn=1000
#PBS -l walltime=168:00:00
#PBS -l mem=8400gb

# Startup scripts

source $HOME/anaconda3/bin/activate
conda activate sc

python $HOME/github/SRAscraper/snakemake_wrapper/scripts/metadata.py

python $HOME/github/SRAscraper/snakemake_wrapper/scripts/fastq_download.py

export PATH=$HOME/cellranger-8.0.1/bin:$PATH

python $HOME/github/SRAscraper/snakemake_wrapper/scripts/cellranger.py

conda activate scpathoquant

python $HOME/github/SRAscraper/snakemake_wrapper/scripts/scpathoquant.py

exit

