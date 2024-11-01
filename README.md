`SRAscraper` (Pipeline for downloading datasets from SRA and GEO databases) First time using `SRAscraper`? We recommend having a look at our [step-by-step guide](https://github.com/MarWoes/wg-blimp/wiki/Tutorial).

#*** WORKING

## Requirements
To run `SRAscraper` you need a UNIX environment that contains a [Bioconda](http://bioconda.github.io/) setup.

## Installation

### Conda environment
First set up a conda environment that will install all the required software needed to run the SRAscraper pipeline. The conda environment can be made from the environment.yaml file included in the github repo after it has been cloned to your local machine.
```
conda env create -f environment.yaml

conda activate SRAscraper
```

###  Install SRAscraper binary
You can install `SRAscraper` from the setup.py into the conda environment by now running the following command in the cloned github repository.
```
pip install .
```

## Running SRAscraper

### SRAscraper pipeline

`SRAscraper` is a cli wrapper for the SRAscraper pipeline implemented using [Snakemake](http://snakemake.readthedocs.io/). In general, a pipeline config is fed to the Snakemake workflow and the corresponding tools are called. However, `wg-blimp` also provides some commands to ease creation of config files, or working without config files altogether.

The command `SRAscraper create-config . config.yaml` is the first step in the pipeline and will generate the `config.yaml` containing all the parameters which will be used to run the pipeline and can be modified easily.

********
Start here
*******

The folder structure created by `wg-blimp run-snakemake` will look as follows:

* alignment - contains all bam/bai files
* dmr - contains dmr files by different callers
* logs - each pipeline step deposits its logs here
* methylation - methylation bedgraph files
* qc - multiqc and other qc related files
* raw - text files describing which fastq files have been used for each sample
* segmentation - methylome segments (UMRs/LMRs/PMDs) as computed by MethylSeekR
* config.yaml - configuration file used for the analysis

It is recommended to check the *raw* folder if all samples contain the correct raw fastq source files.
When in doubt, `wg-blimp` also allows for explicit association of samples and read files by setting `sample_fastq_csv` in the configuration file.
An example csv file could look as follows (column names must be set to `sample`, `forward` and `reverse`):
```
sample,forward,reverse
sample1,/my/path/sample1_L1_1.fq.gz,/my/path/sample1_L1_2.fq.gz
sample1,/my/path/sample1_L2_1.fq.gz,/my/path/sample1_L2_2.fq.gz
sample2,/my/path/sample2_L1_1.fq.gz,/my/path/sample2_L1_2.fq.gz
sample3,/my/path/sample3_L1_1.fq.gz,/my/path/sample3_L1_2.fq.gz
```

## Example

Some example `.fastq` can be found on [Sciebo](https://uni-muenster.sciebo.de/s/7vpqRSEATYcvlnP). You can use the command
```
wg-blimp run-snakemake fastq/ chr22.fasta blood1,blood2 sperm1,sperm2 results --cores=8 --aligner=gembs
```

Please note that the pipeline commands also allow a `--use-sample-files` option so sample groups can be loaded from text files instead of comma separates files.


## Config parameters

The following entries are used for running the Snakemake pipeline and may be specified in the `config.yaml` files:

| Key | Value |
| --- | ----- |
| *aligner* | Aligner to be used by pipeline. Choose either gemBS or bwa-meth. |
| *annotation_allowed_biotypes* | Only genes with this biotype will be annotated in the DMR table (see https://www.gencodegenes.org/pages/biotypes.html ). |
| *annotation_min_mapq* | When annotating coverage, only use reads with a minimum mapping quality |
| *bsseq_local_correct* | Use local correction for bsseq DMR calling. Usually, setting this to FALSE will increase the number of calls. |
| *cgi_annotation_file* | Gzipped csv file used for cg island annotation. Mandatory for MethylSeekR segmentation. Usually downloaded from UCSC Table Browser. |
| *computing_threads* | Number of processors a single job is allowed to use. Remember to use `--cores` parameter for Snakemake. |
| *dmr_tools* | Tools to use for DMR calling. Available: `bsseq`, `camel`, `metilene`
| *group1* | Samples in first group for DMR analysis |
| *group2* | Samples in second group for DMR analysis |
| *gtf_annotation_file* | GTF file used for annotation of genes and promoters. |
| *io_threads* | IO intensive tools virtually reserve this many cores (while actually using only one) to reduce file system IO load. |
| *java_memory_gb* | Gigabytes of RAM to allocate for Java-based tools. If samples are too large, this must be increased to prevent crashes. |
| *methylation_rate_on_chromosomes* | Compute methylation rates for these chromosome during QC |
| *methylseekr_fdr_cutoff* | FDR cutoff for MethylSeekR segmentation. |
| *methylseekr_methylation_cutoff* | Methylation cutoff for MethylSeekR segmentation. |
| *methylseekr_pmd_chromosome* | Chromosome to compute MethylSeekR alpha values for. |
| *min_cov* | Minimum average coverage for methylation calling |
| *min_cpg* | Minimum number of CpGs in a DMR to be called |
| *min_diff* | Minimum average difference between the two groups for DMR calling |
| *output_dir* | Directory containing all files created by the pipeline |
| *promoter_tss_distances* | Distance interval around TSS's to be recognized as promoters in DMR annotation. |
| *rawdir* | Directory containing .fastq files |
| *rawsuffixregex* | The regular expressions to match for paired reads. By default, Illumina naming conventions are accepted. |
| *ref* | .fasta reference file. "Bisulfited" references and BWA indices will be created automatically by bwa-meth) |
| *repeat_masker_annotation_file* | File containing repetitive regions. Usually generated by RepeatMasker and downloaded from UCSC Table Browser. |
| *sample_fastq_csv* | Optional CSV file containing association between samples and read files. The CSV must contain a header with column names `sample`, `forward` and `reverse`. When this option is set, parameters *rawdir* and *rawsuffixregex* are ignored. |
| *samples* | All samples (usually concatenation of group1 and group2) |
| *target_files* | Files to be generated by the Snakemake workflow |
| *temp_dir* | Directory for temporary files. This option may be used for instances where computation node disk space is limited. |

## Reporting errors / Requesting features
If anything goes wrong using `SRAscraper` or any features are missing, feel free to open an issue on this github repo.

## Citing
Please make sure to cite the [github repository](https://github.com/Diako-Lab/SRAscraper) when using SRAscraper for research purposes.
