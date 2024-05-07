"""Collection of functions for setting up the environment for the pipelines.
"""

__author__ = 'Carolin Brune'

################################################################################
# requirements
################################################################################

import click
import os
import requests
import yaml

from datetime import date
from importlib.resources import files
from pathlib import Path
from tqdm import tqdm
from typing import Literal

from refinegems.utility.set_up import download_config as rg_config
from refinegems.utility.io import load_a_table_from_database
from refinegems.utility.databases import initialise_database

################################################################################
# variables
################################################################################

# config keys
# -----------
HQTB_CONFIG_PATH_OPTIONAL = ['media_gap', 'ncbi_map', 'ncbi_dat','biocyc','universal','pan-core']
HQTB_CONFIG_PATH_REQUIRED = ['annotated_genome','full_sequence','model','diamond',
                             'mnx_chem_prop', 'mnx_chem_xref','mnx_reac_prop','mnx_reac_xref',
                             'media_analysis']


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
def download_mnx(dir:str='MetaNetX', chunk_size:int=1024):
    """Download the data needed from the MetaNetX database.

    Args:
        - dir (str, optional): 
            Directory to write the downloaded files to. 
            Defaults to 'MetaNetX'.
        - chunk_size (int, optional): 
            Size of the chunk of data that is loaded into memory during download. 
            Defaults to 1024.
    """

    for mnx_name,mnx_url in MNX_URL_DICT.items():
        r = requests.get(mnx_url, stream=True)
        total = int(r.headers.get('content-length', 0))
        with open(Path(dir,mnx_name), 'wb') as f, tqdm(desc=Path(dir,mnx_name), total=total,
            unit='iB', unit_scale=True, unit_divisor=1024,) as bar:
            for data in r.iter_content(chunk_size=chunk_size):
                size = f.write(data)
                bar.update(size)


# @TODO
#    add part for the cmpb pipeline
def build_data_directories(pipeline: Literal['hqtb','high-quality template based', 
                                         'cmpb', 'carveme modelpolisher based'], 
                           dir:str, chunk_size:int=2048):
    """Set up the necessary directory structure and download files if possible
    for the given pipeline.

    Args:
        - pipeline (Literal['hqtb','high'): 
            For which pipeline the structure should be.
        - dir (str): 
            Parent directory/ Path to write the structure to.
        - chunk_size (int, optional): 
            Chunk size for the download. Defaults to 2048.

    Raises:
        ValueError: Unknown input for parameter pipeline
    """

    match pipeline:
        case 'hqtb' | 'high-quality template based':
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
        case 'cmpb' | 'carveme modelpolisher based':
            # @TODO
            pass
        case _:
            message = f'Unknown input for parameter pipeline: {pipeline}'
            raise ValueError(message)


# ---------------------
# handling config files
# ---------------------

def download_config(filename:str='my_basic_config.yaml', type:Literal['hqtb-basic','hqtb-advanced','hqtb-defaults','media','cmpb']='hqtb basic'):
    """Load a configuration file from the package and save 
    a copy for the user to edit.

    The media config and the config for the cmpb / CarveMe + Modelpolisher based pipeline
    can be downloaded using 'media' and 'cmpb' respectively

    For the htqb / high-quality template based pipeline:

        Depending on the knowledge of the user, either a 'hqtb-basic' or an 'hqtb-advanced' type
        of configuration file can be downloaded (or 'hqtb-defaults' for developers).

    Args:
        - filename (str, optional): 
            Filename/filepath to save the downloaded config file under. 
            Defaults to 'my_basic_config.yaml'.
        - type (Literal['htqb-basic','htqb-advanced','htqb-defaults','media','cmpb'], optional): 
            The type of file to download. 
            Can be 'hqtb-basic', 'hqtb-advanced' or 'hqtb-defaults' or 'media' or 'cmpb'. 
            Defaults to 'hqtb basic'.

    Raises:
        ValueError: Unknown type of config file detected.
    """

    # copy an examplary version of the config file for the user to edit it
    match type:
        # the 'beginner' version
        case 'hqtb-basic':
            config_file = files('specimen.data.config').joinpath('htqb_basic_config_expl.yaml')
            with open(config_file, "r") as cfg_file, open(filename, 'w') as cfg_out:
                for line in cfg_file:
                    cfg_out.write(line)
        # for advanced users
        case 'hqtb-advanced':
            config_file = files('specimen.data.config').joinpath('htqb_advanced_config_expl.yaml')
            with open(config_file, "r") as cfg_file, open(filename, 'w') as cfg_out:
                for line in cfg_file:
                    cfg_out.write(line)
        # for developer: the config with all internal defaults
        case 'hqtb-defaults':
            config_file = files('specimen.data.config').joinpath('htqb_config_default.yaml')
            with open(config_file, "r") as cfg_file, open(filename, 'w') as cfg_out:
                for line in cfg_file:
                    cfg_out.write(line)
        # media config from refinegems
        case 'media':
            rg_config(filename, type='media')
        # for the cmpb pipeline
        case 'cmpb':
            config_file = files('specimen.data.config').joinpath('cmpb_config.yaml')
            with open(config_file, "r") as cfg_file, open(filename, 'w') as cfg_out:
                for line in cfg_file:
                    cfg_out.write(line)
        # type not found
        case _:
            raise ValueError(F'Unknown type of config file detected: {type}')


# hqtb
# ----

# @TODO -> @TEST changes
def validate_hqtb_config(userc:str) -> dict:
    """Validate a user hqtb config file for use in the pipeline.

    Note: currently not everything is checked, mainly the needed files are.

    Args:
        - userc (str): 
            Path to the user configuration file.

    Raises:
        FileNotFoundError: Directory set for config:data:data:direc does not exist.

    Returns:
        dict: 
            The validated, read-in configuration file, nested (read-in yaml file).
    """
    def dict_recursive_combine(dictA:dict, dictB:dict) -> dict:
        """Helper-function for :py:func:`~specimen.util.set_up.validate_hqtb_config` to combine two configuration file.

        Args:
            - dictA (dict): 
                Information from one config file in dict format.
            - dictB (dict): 
                Information from the other config file in dict format.

        Returns:
            dict: 
                The combined information.
        """
        
        if not isinstance(dictB,dict):
            return dictB
        for key in dictA.keys():
            if key in dictB.keys():
                dictA[key] = dict_recursive_combine(dictA[key], dictB[key])
        return dictA


    # @TODO: extent - include more checks
    def dict_recursive_check(dictA:dict, key:str=None):
        """Helper-function for :py:func:`~specimen.util.set_up.validate_hqtb_config` 
        to check if a configuration is valid to run the hight-quality template based pipeline.

        Args:
            - dictA (dict): 
                Current dictionary or value to be validated.
            - key (str, optional): 
                key of dictA, if it was an entry of a dictionary. 
                Defaults to None.

        Raises:
            TypeError: Missing a required argument in the config file.

            FileNotFoundError: Path does not exist: {dictA}

            FileNotFoundError: Path does not exist: {dictA}
        """

        if not isinstance(dictA,dict):
            if dictA == '__USER__':
                raise TypeError(F'Missing a required argument in the config file: {key}')
            elif key in HQTB_CONFIG_PATH_REQUIRED:
                if os.path.isfile(dictA):
                    return
                else:
                    raise FileNotFoundError(F'Path does not exist: {dictA}')
            elif key in HQTB_CONFIG_PATH_OPTIONAL:
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


    # validate a user config file by checking for missing input
    # by combining it with a default config

    # load both files
    defaultc_path = files('specimen.data.config').joinpath('htqb_config_default.yaml')
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


# cmpb
# ----

# @TODO
#    is it still correct after all the changes???
def save_cmpb_user_input(configpath: str) -> dict[str: str]:
    """This aims to collect user input from the command line to create a cmpb config file, 
    will also save the user input to a config if no config was given

    Args:
        - configpath (str): 
            Path to config file if present
        
    Returns:
        dict: 
            Either loaded config file or created from user input
    """
    if os.path.isfile(configpath):
        with open(configpath) as f:
            config = yaml.safe_load(f)
        print(config)
        return config
    else:
        print('No config or no valid config given, you will be asked for input')
        user_input = {}
        
        update_db = click.confirm('Do you want to update the database?')
        user_input['db_update'] = update_db
        
        if update_db:
            initialise_database()

        out_path = click.confirm('Do you want to keep the output path "../rg_out/"?', default=True)
        if not out_path:
            user_input['out_path'] = click.prompt('Enter your desired output path')
        else:
            user_input['out_path'] = '../rg_out/'
        
        user_input['visualize'] = click.confirm('Do you want to generate visualizations of your model(s)?')
            
        growth_basis = click.prompt('Enter the base uptakes for growth simulation (d for default_uptake, m for minimal_uptake)')
        if growth_basis == 'd':
            user_input['growth_basis'] = 'default_uptake'
        if growth_basis == 'm':
            user_input['growth_basis'] = 'minimal_uptake'
            
        user_input['anaerobic_growth'] = click.confirm('Do you want to simulate anaerobic growth?')
        
        multiple = click.confirm('Do you want to simulate and compare multiple models?')
        user_input['multiple'] = multiple
        if multiple:
            list_of_models = []
            while True:
                file_path = click.prompt('Enter file path to model (or "stop" to stop)')
                if file_path.lower() == 'stop':
                    break
                elif os.path.isfile(file_path):
                    list_of_models.append(file_path)
                    print('Added file:', file_path)
                else:
                    print('File does not exist. Please enter a valid file path.')
            print('The following models will be compared:')
            print(list_of_models)
            user_input['multiple_paths'] = list_of_models
        possible_media = load_a_table_from_database('medium', False)['name'].to_list()
        possible_media_str = '|'.join(possible_media)
        list_of_media = []
        while True:
            medium = click.prompt(f'Enter medium to simulate growth on ({possible_media_str}) (or "stop" to stop)')
            if medium.lower() == 'stop':
                break
            elif medium in possible_media:
                if medium not in list_of_media:
                    list_of_media.append(medium)
                else:
                    print(medium + ' is already in the list.')
            else:
                print('Please choose a medium from the given list.')
        user_input['media'] = list_of_media
        
        single = click.confirm('Do you want to investigate or curate a single model?')
        user_input['single'] = single
        if single:
            not_valid = True
            while not_valid:
                model = click.prompt('Path to your model file.')
                if os.path.isfile(model):
                    user_input['model'] = model
                    not_valid = False
                else:
                    print('File does not exist. Please enter a valid file path')
            user_input['memote'] = click.confirm('Do you want to run MEMOTE (takes some time)?')    
            user_input['modelseed'] = click.confirm('Do you want to compare your model entities to the ModelSEED database?')
        
            gap_analysis = click.confirm('Do you want to run a gap analysis?') 
            user_input['gap_analysis'] = gap_analysis
            if gap_analysis:
                gap_analysis_params = {}
                db_to_compare = click.prompt('One of the choices KEGG|BioCyc|KEGG+BioCyc') #|GFF
                gap_analysis_params['db_to_compare'] = db_to_compare
                if db_to_compare == 'KEGG' or db_to_compare == 'KEGG+BioCyc':
                    user_input['organismid'] = click.prompt('Enter the KEGG Organism code')
                    user_input['gff_file'] = click.prompt('Enter the path to your organisms RefSeq GFF file')
                if db_to_compare == 'BioCyc' or db_to_compare == 'KEGG+BioCyc':
                    Path0 = click.prompt('Enter the path to your BioCyc TXT file containing a SmartTable with the columns \'Accession-2\' and \'Reaction of gene\'')
                    Path1 = click.prompt('Enter the path to your BioCyc TXT file containing a SmartTable with all reaction relevant information')
                    Path2 = click.prompt('Enter the path to your Biocyc TXT file containing a SmartTable with all metabolite relevant information')
                    Path3 = click.prompt('Enter path to protein FASTA file used as input for CarveMe')
                    gap_analysis_params['biocyc_files'] = [Path0, Path1, Path2, Path3]
                user_input['gap_analysis_params'] = gap_analysis_params
                
            mod = click.confirm('Do you want to use functions to modify your model?')
            if mod:
                
                new_path = click.confirm('Do you want to save your modified model to ' + user_input['out_path'] + '<model.id>_modified_<today>.xml?')
                if new_path:
                    user_input['model_out'] = 'stdout'
                else:
                    user_input['model_out'] = click.prompt('Enter path and filename to where to save the modified model')
                
                gapfill_model = click.confirm('Do you want to fill gaps in your model?')
                user_input['gapfill_model'] = gapfill_model
                
                if gapfill_model:
                    if not gap_analysis:
                        user_input['gap_analysis_file'] = click.prompt('Enter path to Excel file with which gaps should be filled')
                
                user_input['keggpathways'] = click.confirm('Do you want to add KEGG Pathways?')
                    
                user_input['sboterms'] = click.confirm('Do you want to update the SBO Terms?')
                
                user_input['charge_corr'] = click.confirm('Do you want to add charges to uncharged metabolites?')
                    
                man_cur = click.confirm('Do you want to modify your model with the manual curations table?')
                user_input['man_cur'] = man_cur

                if man_cur:
                    entrez_email = click.prompt('Email to access NCBI Entrez')
                    user_input['entrez_email'] = entrez_email
                    man_cur_type = click.prompt('Enter type of curation (gapfill|metabs)')
                    user_input['man_cur_type'] = man_cur_type
                    man_cur_table = click.prompt('Enter the path to the manual curations table')
                    user_input['man_cur_table'] = man_cur_table

                polish = click.confirm('Do you want to polish the model?')
                user_input['polish'] = polish

                if polish:
                    entrez_email = click.prompt('Email to access NCBI Entrez')
                    user_input['entrez_email'] = entrez_email
                    id_db = click.prompt('What database is your model based on? BIGG|VMH')
                    user_input['id_db'] = id_db
                    gff_file = click.prompt('If possible, provide the path to the RefSeq GFF file of your organism')
                    user_input['gff_file'] = gff_file if gff_file != 'None' else None
                    protein_fasta = click.prompt('If possible, provide the path to the Protein FASTA file used for CarveMe')
                    user_input['protein_fasta'] = protein_fasta if protein_fasta != 'None' else None
                    lab_strain = not click.confirm('Does your modeled organism have a database entry?', default=True)
                    user_input['lab_strain'] = lab_strain
                    organismid = click.prompt('If possible, provide the KEGG organism code of your organism')
                    user_input['organismid'] = organismid if (organismid != 'None') else None
                    
                biomass = click.confirm('Do you want to check & normalise the biomass function(s)?')
                user_input['biomass'] = biomass
                    
            else:
                user_input['keggpathways'] = False
                user_input['polish'] = False
                user_input['biomass'] = False
                user_input['sboterms'] = False
                user_input['charge_corr'] = False
                user_input['gapfill_model'] = False
                user_input['man_cur'] = False
            
        today = date.today().strftime("%Y%m%d")
        
        print('This is your input:')
        print(user_input)
        if not os.path.isdir(user_input['out_path']):
            print('Given out_path is not yet a directory, creating ' + user_input['out_path'])
            os.makedirs(user_input['out_path'])
        with open(user_input['out_path'] + 'user_input_' + str(today) + '.yaml', 'w') as f:
            yaml.dump(user_input, f)
        print('Your input was saved as yaml to '+ user_input['out_path'] + 'user_input_' + str(today) + '.yaml')
        return user_input

