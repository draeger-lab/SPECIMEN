"""Collection of functions for setting up the environment for the pipeline.
"""

__author__ = 'Carolin Brune'

################################################################################
# requirements
################################################################################

from importlib.resources import files
import os
from pathlib import Path
import requests
from tqdm import tqdm
import yaml

################################################################################
# variables
################################################################################

# config keys
# -----------
CONFIG_PATH_OPTIONAL = ['media', 'ncbi_map', 'ncbi_dat','biocyc','universal','pan-core']
CONFIG_PATH_REQUIRED = ['annotated_genome','full_sequence','model','diamond',
                        'mnx_chem_prop', 'mnx_chem_xref','mnx_reac_prop','mnx_reac_xref']


# external databases
# ------------------
MNX_CHEM_XREF_URL = 'https://www.metanetx.org/cgi-bin/mnxget/mnxref/chem_xref.tsv'
MNX_CHEM_PROP_URL = 'https://www.metanetx.org/cgi-bin/mnxget/mnxref/chem_prop.tsv'
MNX_REAC_XREF_URL = 'https://www.metanetx.org/cgi-bin/mnxget/mnxref/reac_xref.tsv'
MNX_REAC_PROP_URL = 'https://www.metanetx.org/cgi-bin/mnxget/mnxref/reac_prop.tsv'
MNX_URL_DICT = {'chem_prop.tsv':MNX_CHEM_PROP_URL, 'chem_xref.tsv':MNX_CHEM_XREF_URL,
                'reac_prop.tsv':MNX_REAC_PROP_URL, 'reac_xref.tsv':MNX_REAC_XREF_URL}

################################################################################
# functions
################################################################################

# ----------------------
# setup data (structure)
# ----------------------
# @TEST : deleted BiGG part, since its already covered with refinegems

def download_mnx(dir='MetaNetX/', chunk_size=1024):
    """Download the data needed from the MetaNetX database.

    :param dir: Directory to write the diwnloaded files to.
    :type dir: string, should end with '/'.
    :param chunk_size: Size of the chunk of data that is loaded into memory during download.
        Default is 1024.
    :type chunk_size: int
    """

    for mnx_name,mnx_url in MNX_URL_DICT.items():
        r = requests.get(mnx_url, stream=True)
        total = int(r.headers.get('content-length', 0))
        with open(Path(dir,mnx_name), 'wb') as f, tqdm(desc=Path(dir,mnx_name), total=total,
            unit='iB', unit_scale=True, unit_divisor=1024,) as bar:
            for data in r.iter_content(chunk_size=chunk_size):
                size = f.write(data)
                bar.update(size)


def build_data_directories(dir, chunk_size=2048):
    """Set up the directory structure for the data and download the files
    from MetaNetX.

    :param dir: Parent folder to write the subfolder structure to.
    :type dir: string
    :param chunk_size: Size of the chunks of data while downloading.
        Default is 2048. Only used for downloading the MetaNetX files.
    :type chunk_size: int
    """

    # create the data directory structure
    print('Creating directory structure...')
    DATA_DIRECTORIES = ['annotated_genomes', 'BioCyc', 'RefSeqs',
                        'medium', 'MetaNetX', 'pan-core-models', 'template-models',
                        'universal-models']
    for sub_dir in DATA_DIRECTORIES:
        new_dir = Path(dir,sub_dir)
        try:
            Path(new_dir).mkdir(parents=True, exist_ok=False)
            print(F'Creating new directory {new_dir}')
        except FileExistsError:
            print(F'Directory {new_dir} already exists.')

    # download data for those directories where this is possible
    print('Downloading MetaNetX...')
    download_mnx(Path(dir,'MetaNetX'), chunk_size=chunk_size)


# ---------------------
# handling config files
# ---------------------

def download_config(filename='./my_basic_config.yaml', type='basic'):
    """Load a configuration file from the package and save a copy for the user to edit.

    Depending on the knowledge of the user, either a 'basic' or an 'advanced' type
    of configuration file can be downloaded.

    :param filename: Filename/filepath to save the downloaded config file under.
    :type filename: string
    :param type: The type of file to download. Can either be 'basic' or 'advanced'.
        Default is 'basic'.
    :type type: string

    :raises: :class:`ValueError`: 'Unknown type of config file detected'
    """


    # copy an examplary version of the config file for the user to edit it
    match type:
        # the 'beginner' version
        case 'basic':
            config_file = files('specimen.data.config').joinpath('basic_config_expl.yaml')
            with open(config_file, "r") as cfg_file, open(filename, 'w') as cfg_out:
                for line in cfg_file:
                    cfg_out.write(line)
        # for advanced users
        case 'advanced':
            config_file = files('specimen.data.config').joinpath('advanced_config_expl.yaml')
            with open(config_file, "r") as cfg_file, open(filename, 'w') as cfg_out:
                for line in cfg_file:
                    cfg_out.write(line)
        # type not found
        case _:
            raise ValueError(F'Unknown type of config file detected: {type}')


def dict_recursive_combine(dictA, dictB):
    """Helper-function to combine two configuration file.

    """
    if not isinstance(dictB,dict):
        return dictB
    for key in dictA.keys():
        if key in dictB.keys():
            dictA[key] = dict_recursive_combine(dictA[key], dictB[key])
    return dictA


# @TODO: extent - include more checks
def dict_recursive_check(dictA, key=None):
    """Function to check if a configuration is valid to run the pipeline.

    :param dictA: Current dictionary or value to be validated.
    :type dictA: dict or other
    :param key: key of dictA, if it was an entry of a dictionary.
        Default is None.
    :type key: string
    """

    if not isinstance(dictA,dict):
        if dictA == '__USER__':
            raise TypeError(F'Missing a required argument in the config file: {key}')
        elif key in CONFIG_PATH_REQUIRED:
            if os.path.isfile(dictA):
                return
            else:
                raise FileNotFoundError(F'Path does not exist: {dictA}')
        elif key in CONFIG_PATH_OPTIONAL:
            if dictA and os.path.isfile(dictA):
                return
            elif not dictA:
                return
            else:
                raise FileNotFoundError(F'Path does not exist: {dictA}')
        else:
            pass
            # @TODO
            #   extend function to cover more of the options
        return
    else:
        for key in dictA.keys():
            dict_recursive_check(dictA[key], key)
    return


# @TODO need to fix this after fixing the config
def validate_config(userc):
    """Validate a user config file for use in the pipeline.

    Note: currently not everything is checked, mainly the needed files are.

    :param userc: Path to the user configuration file.
    :type userc: string

    :raise: :FileNotFoundError: Directory set for config:data:data:direc does not exist.

    :returns: The validated, read-in configuration file.
    :rtype: dict, nested (read-in yaml file)
    """

    # validate a user config file by checking for missing input
    # by combining it with a default config

    # load both files
    defaultc_path = files('config').joinpath('config_default.yaml')
    with open(defaultc_path, "r") as cfg_def, open(userc, 'r') as cfg_usr:
        config_d = yaml.load(cfg_def, Loader=yaml.loader.FullLoader)
        config_u = yaml.load(cfg_usr, Loader=yaml.loader.FullLoader)

    # combine
    combined_config = dict_recursive_combine(config_d, config_u)

    # check for missing or problematic values
    if combined_config['data']['data_direc']:
        if os.path.isdir(combined_config['data']['data_direc']):
            for key in combined_config['data']:
                if combined_config['data'][key] and key != 'data_direc':
                    combined_config['data'][key] = combined_config['data']['data_direc']  + combined_config['data'][key]
            dict_recursive_check(combined_config, key=None)
        else:
            raise FileNotFoundError('Directory set for config:data:data_direc does not exist.')
    else:
        dict_recursive_check(combined_config, key=None)

    return combined_config


