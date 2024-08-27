import os
import requests
from ruamel.yaml import YAML

script_dir = os.path.dirname(os.path.realpath(__file__))

DEFAULT_OPTIONALS_FILE = os.path.join(script_dir, 'optionals.yaml')

yaml = YAML(typ='safe')
yaml.default_flow_style = False

def get_default_optional_parameters():

    with open(DEFAULT_OPTIONALS_FILE) as f:

        default_optionals = yaml.load(f)

        return default_optionals

def get_default_config(genome_build, cores_per_job, nodes, output_dir):

    mandatory_parameters = {
        'output_dir': os.path.abspath(output_dir),
        'computing_threads': cores_per_job,
        'nodes': nodes
    }

    optional_parameters = get_default_optional_parameters()

    return {
        **mandatory_parameters,
        **optional_parameters
    }


def dump_config(config_dict, target_file):

    with (open(target_file, 'w')) as f:

        yaml.dump(config_dict, f)


def create_config(genome_build, cores_per_job, nodes, output_dir, target_yaml):

    config_yaml = get_default_config(genome_build, cores_per_job, nodes, output_dir)

    dump_config(config_yaml, target_yaml)
