#!/bin/bash

#PBS -N Pile_of_viruses
#PBS -l nodes=1:ppn=1000
#PBS -l walltime=168:00:00
#PBS -l mem=8400gb

# Startup scripts

source $HOME/anaconda3/bin/activate
conda activate sc

python $HOME/github/SRAscraper/snakemake_wrapper/scripts/remake_virus/remake_viral_ref.py

exit

