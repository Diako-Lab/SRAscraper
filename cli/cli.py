import click
import cli.config
import cli.snakemake


@click.group()
def main():
    pass


@main.command(help='Create a config YAML file for running the OnmiScrape pipeline.')
@click.option('--genome_build', type=click.Choice(['hg19','hg38', 'mmul10', 'None']), default='hg38', help='Build of the reference used for annotation.')
@click.option('--cores-per-job', type=click.INT, default='1', help='The number of cores to use per job.')
@click.option('--nodes', type=click.INT, default=1, help='Number of nodes in a cluster to use in cluster mode. SERVER ONLY.')
@click.argument('output_dir')
@click.argument('target_yaml')
def create_config(genome_build, cores_per_job, nodes, output_dir, target_yaml):

    cli.config.create_config(genome_build, cores_per_job, nodes, output_dir, target_yaml)


@main.command(help='Run the OnmiScrape pipeline using a config file.')
@click.option('--dry-run', is_flag=True, default=False, help='Only dry-run the workflow.')
@click.option('--cores', required=True, type=click.INT, help='The maximum number of cores to use for running the pipeline. Cores per job are set in configuration file. In cluster mode, this sets cores used per node.')
@click.option('--cluster', default=None, type=click.STRING, help='Submission command snakemake uses for cluster usage. Setting this parameter enables snakemake\'s cluster mode.')
@click.option('--nodes', default=1, type=click.INT, help='Number of nodes to use in cluster mode. SERVER ONLY')
@click.argument('config_yaml')
def run_snakemake_from_config(dry_run, cores, cluster, nodes, config_yaml):

    cli.OmniScrape.run_from_config(dry_run, cores, cluster, nodes, config_yaml)

