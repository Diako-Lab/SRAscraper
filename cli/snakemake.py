import os
import snakemake
import yaml

# get right directory, see https://stackoverflow.com/questions/4934806/how-can-i-find-scripts-directory-with-python
script_dir = os.path.dirname(os.path.realpath(__file__))
snakefile_location = os.path.join(script_dir,'..','snakemake_wrapper', 'Snakefile')
conda_prefix = os.path.abspath(os.path.join(script_dir, '..', 'snakemake_wrapper', 'conda'))


def run_config(dry_run, config_yaml):

    print("[INFO] Invoking Snakemake with config {}".format(config_yaml))

    with open(config_yaml, 'r') as stream:
        data_loaded = yaml.safe_load(stream)

    output_dir = data_loaded['output_dir']

    cores = data_loaded['computing_threads']
    
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)
    
    finished_successfully = snakemake.snakemake(
        snakefile=snakefile_location,
        configfiles=[config_yaml],
        cores=cores,
        dryrun=dry_run,
        printshellcmds=True,
        use_conda=True,
        conda_prefix=conda_prefix
    )

    if not finished_successfully:
        os.sys.exit(os.EX_SOFTWARE)


