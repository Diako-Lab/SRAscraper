import os, io
import requests
import yaml

script_dir = os.path.dirname(os.path.realpath(__file__))

DEFAULT_OPTIONALS_FILE = os.path.join(script_dir, 'optionals.yaml')

def check_file_path(file_path):
    """Checks if a file path is valid and returns an error if not."""

    if not os.path.exists(os.path.abspath(os.path.realpath(file_path))):
        raise FileNotFoundError(f"File not found: {file_path}")

    return os.path.abspath(os.path.realpath(file_path)) 

def get_default_optional_parameters():

    with open(DEFAULT_OPTIONALS_FILE, 'r') as f:

        default_optionals = yaml.safe_load(f)

        return default_optionals

def get_default_config(cores, output_dir, NCBI_search_txt):

    mandatory_parameters = {
        'output_dir': os.path.abspath(output_dir),
        'NCBI_search_txt': check_file_path(NCBI_search_txt),
        'computing_threads': cores
    }

    optional_parameters = get_default_optional_parameters()

    return {
        **mandatory_parameters,
        **optional_parameters
    }


def dump_config(config_dict, target_file):

    with (io.open(target_file, 'w', encoding='utf8')) as f:

        yaml.dump(config_dict, f, default_flow_style=False, allow_unicode=True)


def create_config(cores, output_dir, NCBI_search_txt, target_yaml):

    config_yaml = get_default_config(cores, output_dir, NCBI_search_txt)

    dump_config(config_yaml, target_yaml)


