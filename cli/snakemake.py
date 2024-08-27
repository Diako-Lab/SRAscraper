import os
import snakemake
import cli.config

# get right directory, see https://stackoverflow.com/questions/4934806/how-can-i-find-scripts-directory-with-python
script_dir = os.path.dirname(os.path.realpath(__file__))
snakefile_location = os.path.join(script_dir,'..','snakemake_wrapper', 'Snakefile')
conda_prefix = os.path.abspath(os.path.join(script_dir, '..', 'snakemake_wrapper', 'conda'))

def run_from_config(dry_run, config_yaml, cores, cluster_command, nodes):

    print("[INFO] Invoking Snakemake with config {}, {} cores, and {} nodes.".format(config_yaml, cores, nodes))

    finished_successfully = snakemake.snakemake(
        snakefile=snakefile_location,
        configfiles=[config_yaml],
        dryrun=dry_run,
        cores=cores,
        local_cores=cores,
        nodes=nodes,
        printshellcmds=True,
        use_conda=True,
        conda_prefix=conda_prefix,
        cluster=cluster_command
    )
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)    


    if not finished_successfully:
        os.sys.exit(os.EX_SOFTWARE)


def create_cofig(genome_build, cores, nodes, output_dir, yaml_config_file):

    if not os.path.exists(output_dir):
        os.mkdir(output_dir)


        yaml_config_file = os.path.join(output_dir, 'config.yaml')
