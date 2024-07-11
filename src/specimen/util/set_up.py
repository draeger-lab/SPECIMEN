"""Collection of functions for setting up the environment for the pipelines.
"""

__author__ = 'Carolin Brune'

################################################################################
# requirements
################################################################################

import click
import logging
import os
import requests
import yaml

from datetime import date
from importlib.resources import files
from pathlib import Path
from tqdm import tqdm
from typing import Literal,Union

from refinegems.utility.set_up import download_config as rg_config

################################################################################
# variables
################################################################################

# config keys
# -----------

# config keys for pipeline files
HQTB_CONFIG_PATH_OPTIONAL = ['media_gap', 'ncbi_map', 'ncbi_dat','biocyc','universal','pan-core'] #: :meta: 
HQTB_CONFIG_PATH_REQUIRED = ['annotated_genome','full_sequence','model','diamond',
                             'mnx_chem_prop', 'mnx_chem_xref','mnx_reac_prop','mnx_reac_xref',
                             'media_analysis'] #: :meta: 
CMPB_CONFIG_PATHS_REQUIRED = ['annotated_genome','mediapath','modelpath','refseq_gff','protein_fasta',
                     'biocyc_files','full_genome_sequence'] #: :meta: 
CMPB_CONFIG_PATHS_OPTIONAL = [] # :meta:
PIPELINE_PATHS_OPTIONAL = {'hqtb':HQTB_CONFIG_PATH_OPTIONAL,
                  'cmpb':CMPB_CONFIG_PATHS_OPTIONAL} #: :meta: 
PIPELINE_PATHS_REQUIRED = {'hqtb':HQTB_CONFIG_PATH_REQUIRED,
                  'cmpb':CMPB_CONFIG_PATHS_REQUIRED} #: :meta: 
# config keys for pipelines directories
PIPELINE_DIR_PATHS = ['dir']


# external databases
# ------------------
MNX_CHEM_XREF_URL = 'https://www.metanetx.org/cgi-bin/mnxget/mnxref/chem_xref.tsv' #: :meta: 
MNX_CHEM_PROP_URL = 'https://www.metanetx.org/cgi-bin/mnxget/mnxref/chem_prop.tsv' #: :meta: 
MNX_REAC_XREF_URL = 'https://www.metanetx.org/cgi-bin/mnxget/mnxref/reac_xref.tsv' #: :meta: 
MNX_REAC_PROP_URL = 'https://www.metanetx.org/cgi-bin/mnxget/mnxref/reac_prop.tsv' #: :meta: 
MNX_URL_DICT = {'chem_prop.tsv':MNX_CHEM_PROP_URL, 'chem_xref.tsv':MNX_CHEM_XREF_URL,
                'reac_prop.tsv':MNX_REAC_PROP_URL, 'reac_xref.tsv':MNX_REAC_XREF_URL} #: :meta: 

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
        - ValueError: Unknown input for parameter pipeline
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

    For the hqtb / high-quality template based pipeline:

        Depending on the knowledge of the user, either a 'hqtb-basic' or an 'hqtb-advanced' type
        of configuration file can be downloaded (or 'hqtb-defaults' for developers).

    Args:
        - filename (str, optional): 
            Filename/filepath to save the downloaded config file under. 
            Defaults to 'my_basic_config.yaml'.
        - type (Literal['hqtb-basic','hqtb-advanced','hqtb-defaults','media','cmpb'], optional): 
            The type of file to download. 
            Can be 'hqtb-basic', 'hqtb-advanced' or 'hqtb-defaults' or 'media' or 'cmpb'. 
            Defaults to 'hqtb basic'.

    Raises:
        - ValueError: Unknown type of config file detected.
    """

    # copy an examplary version of the config file for the user to edit it
    match type:
        # the 'beginner' version
        case 'hqtb-basic':
            config_file = files('specimen.data.config').joinpath('hqtb_basic_config_expl.yaml')
            with open(config_file, "r") as cfg_file, open(filename, 'w') as cfg_out:
                for line in cfg_file:
                    cfg_out.write(line)
        # for advanced users
        case 'hqtb-advanced':
            config_file = files('specimen.data.config').joinpath('hqtb_advanced_config_expl.yaml')
            with open(config_file, "r") as cfg_file, open(filename, 'w') as cfg_out:
                for line in cfg_file:
                    cfg_out.write(line)
        # for developer: the config with all internal defaults
        case 'hqtb-defaults':
            config_file = files('specimen.data.config').joinpath('hqtb_config_default.yaml')
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
def validate_config(userc:str, pipeline:Literal['hqtb','cmpb']='hqtb') -> dict:
    """Validate a user hqtb config file for use in the pipeline.

    .. note::
    
        Currently not everything is checked, mainly the needed files are.

    .. warning::

        Currently working not perfectly, under construction

    Args:
        - userc (str): 
            Path to the user configuration file.

    Raises:
        - FileNotFoundError: Directory set for config:data:data:direc does not exist.

    Returns:
        dict: 
            The validated, read-in configuration file, nested (read-in yaml file).
    """
    def dict_recursive_combine(dictA:dict, dictB:dict) -> dict:
        """Helper-function for :py:func:`~specimen.util.set_up.validate_config` to combine two configuration file.

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
    
    def dict_recursive_overwrite(dictA:dict) -> dict:
        """Helper-function for :py:func:`~specimen.util.set_up.validate_config` to combine two configuration file.

        Args:
            - dictA (dict): 
                The dictionary to validate

        Raises:
            - TypeError: Missing file/path

        Returns:
            dict: 
                the Dictionary with USER overwritten as None
        """

        if not isinstance(dictA,dict):
            # check for missing input
            if dictA == '__USER__':
                raise TypeError(F'Missing a required argument in the config file')
            elif dictA == 'USER':
                # @TODO 
                mes = 'Keyword USER detected in config. Either due to skipped options or missing required information.\nReminder: this may lead to downstream problems.'
                logging.warning(mes)
                return None
            else:
                return dictA

        for key in dictA.keys():
            dictA[key] = dict_recursive_overwrite(dictA[key])
        return dictA


    # @TODO: extent - include more checks
    # @TEST: does the current version work correctly? (or even work at all) 
    def dict_recursive_check(dictA:dict, key:str=None, 
                             pipeline:Literal['hqtb','cmpb']='hqtb'):
        """Helper-function for :py:func:`~specimen.util.set_up.validate_config` 
        to check if a configuration is valid to run the hight-quality template based pipeline.

        Args:
            - dictA (dict): 
                Current dictionary or value to be validated.
            - key (str, optional): 
                key of dictA, if it was an entry of a dictionary. 
                Defaults to None.

        Raises:
            - TypeError: Missing a required argument in the config file.
            - FileNotFoundError: Path does not exist: {dictA}
            - FileNotFoundError: Path does not exist: {dictA}
        """

        if not isinstance(dictA,dict):
            # required file paths
            if key in PIPELINE_PATHS_REQUIRED[pipeline]: 
                if dictA and os.path.isfile(dictA):
                    return
                else: 
                    raise FileNotFoundError(F'Path does not exist: {dictA}')
            # missing optional file
            elif key in PIPELINE_PATHS_OPTIONAL and not dictA:
                # @TODO, already warning in Overwrite?
                return 
            # optional file paths
            elif key in PIPELINE_PATHS_OPTIONAL:
                if isinstance(str,dictA):
                    if os.path.isfile(dictA):
                        return
                    else: 
                        raise FileNotFoundError(F'Path does not exist: {dictA}')
                if isinstance(list,dictA):
                    for entry in dictA:
                        if dictA and os.path.isfile(dictA):
                            return
                        elif not dictA:
                            # @TODO
                            pass
                        else: 
                            raise FileNotFoundError(F'Path does not exist: {dictA}')
            elif key in PIPELINE_DIR_PATHS:
                if dictA and os.path.exists(dictA):
                    return
                else:
                    raise FileNotFoundError(F'Directory does not exist: {dictA}')
            # not found or missing
            else:
                # @TODO
                pass
                
            return

        else:
            for key in dictA.keys():
                dict_recursive_check(dictA[key], key, pipeline)
        return


    # validate a user config file by checking for missing input
    # by combining it with a default config

    # load both files
    match pipeline:
        case 'hqtb':
            defaultc_path = files('specimen.data.config').joinpath('hqtb_config_default.yaml')
        case 'cmpb':
            defaultc_path = files('specimen.data.config').joinpath('cmpb_config.yaml')
        case _:
            raise ValueError(f'Unknown input for pipeline: {pipeline}')
    
    with open(defaultc_path, "r") as cfg_def, open(userc, 'r') as cfg_usr:
        config_d = yaml.load(cfg_def, Loader=yaml.loader.FullLoader)
        config_u = yaml.load(cfg_usr, Loader=yaml.loader.FullLoader)

    # combine
    combined_config = dict_recursive_combine(config_d, config_u)

    # overwrite __USER__ and USER
    combined_config = dict_recursive_overwrite(combined_config)

    # check for missing or problematic values
    # special case for HQTB pipeline with relative paths
    if 'data' in combined_config.keys() and 'data_direc' in combined_config['data'].keys() and combined_config['data']['data_direc']:
        if os.path.isdir(combined_config['data']['data_direc']):
            for key in combined_config['data']:
                if combined_config['data'][key] and key != 'data_direc':
                    combined_config['data'][key] = combined_config['data']['data_direc']  + combined_config['data'][key]
            dict_recursive_check(combined_config, key=None, pipeline=pipeline)
        else:
            raise FileNotFoundError('Directory set for config:data:data_direc does not exist.')
    # normal recursion for validation
    else:
        dict_recursive_check(combined_config, key=None, pipeline=pipeline)

    return combined_config


# cmpb
# ----

# @TEST
def save_cmpb_user_input(configpath:Union[str,None]=None) -> dict:
    """Guide the user step by step throuh the creation of the configuration for a cmpb pipeline run
    (via commandline).

    Args:
        - configpath (Union[str,None], optional): 
            Path to a file to save the config under. Defaults to None.

    Returns:
        dict: 
            The configuration in dictionary format.
    """

        
    print('No config or no valid config given, you will be asked for input')

    config_file = files('specimen.data.config').joinpath('cmpb_config.yaml')
    with open(config_file, "r") as cfg:
        config = yaml.load(cfg, Loader=yaml.loader.FullLoader)

    # if model, get path
    has_model = click.prompt('Do you already have a draft model e.g. created with CarveMe?', type=click.Choice(['y','n']), show_choices=True)
    match has_model:
        case 'y':
            modelpath = click.prompt('Enter the path to your model', type=click.Path(exists=True))
            config['input']['modelpath'] = modelpath
        case 'n':
            pass
     
    # required input
    # --------------
    print('------------')
    print('The following information is REQUIRED for the pipeline')
    print('------------')
    config['input']['annotated_genome'] = click.prompt('Enter the path to your annotated genome file', type=click.Path(exists=True))
    config['input']['mediapath'] = click.prompt('Enter the path to a media configuration file for growth simulation', type=click.Path(exists=True))

    # general options
    # ---------------
    print('------------')
    print('General options')
    print('------------')

    # output directory
    config['general']['dir'] = click.prompt('Enter your desired output directory path', type=click.Path())
    
    # colour 
    set_col = click.prompt('Do you want to use the default colour map YlGn for the visualisation?', type=click.Choice(['y','n']), show_choices=True)
    match set_col:
        case 'n':
            colours = click.prompt('Enter your choosen colour scheme', type=str)
            config['general']['colours'] = colours
        case 'y':
            pass

    # save all models or not
    save_models = click.prompt('Do you want to save the model separatly after each step?', type=click.Choice(['y','n']), show_choices=True)
    match save_models:
        case 'y':
            config['general']['save_all_models'] = True
        case 'n':
            config['general']['save_all_models'] = False

    # run memote always y/n
    run_memote = click.prompt('Do you want run memote after each step?', type=click.Choice(['y','n']), show_choices=True)
    match run_memote:
        case 'y':
            config['general']['memote_always_on'] = True
        case 'n':
            config['general']['memote_always_on'] = False
    
    # run stats always y/n
    models_stats = click.prompt('Do you want run memote after each step?', type=click.Choice(['y','n']), show_choices=True)
    match models_stats:
        case 'y':
            config['general']['stats_always_on'] = True
        case 'n':
            config['general']['stats_always_on'] = False

    # some additional, sometimes required, sometimes optional files
    refseq = click.prompt('If you want to run a gap analysis with KEGG or have a CarveMe model, please enter the path to your refseq gff file', type=click.Path(exists=True))
    config['general']['refseq_organism_id'] = refseq

    kegg_org_id = click.prompt('If you want to run a gap analysis with KEGG or have a CarveMe model, please enter the path to your refseq gff file', type=click.Path(exists=True))
    config['general']['refseq_organism_id'] = kegg_org_id

    # part-specific 
    # -------------
    print('------------')
    print('Part-specific options')
    print('------------')

    # model polish
    carveme = click.prompt('Is your draft model CarveMe-based or will CarveMe be run?', type=click.Choice(['y','n']), show_choices=True)
    if carveme == 'y':
        email = click.prompt('Enter an email address for the connection to NCBI Entrez', type=str)
        config['cm-polish']['email'] = email
        labs = click.prompt('Do you have a lab strain?', type=click.Choice(['y','n']), show_choices=True)
        labs = True if labs == 'y' else False
        config['cm-polish']['is_lab_strain'] = labs
        if labs:
            prot_fa = click.prompt('Please enter the path to your protein fasta file', type=click.Path(exists=True))
            config['cm-polish']['protein_fasta'] = prot_fa
        else:
            check_prot_fa = click.prompt('Do you want to/can provide the protein fasta?', type=click.Choice(['y','n']), show_choices=True)
            if check_prot_fa == 'y':
                prot_fa = click.prompt('Please enter the path to your protein fasta file', type=click.Path(exists=True))
                config['cm-polish']['protein_fasta'] = prot_fa
            else:
                config['cm-polish']['protein_fasta'] = None
    

    # gapfilling
    gap_analysis = click.prompt('Do you want to run a gap analysis?', type=click.Choice(['y','n']), show_choices=True) 

    if gap_analysis == 'y':
        db_to_compare = click.prompt('Available database information for the gapfilling', type=click.Choice(['KEGG','BioCyc','KEGG+BioCyc']), show_choices=True) 
        config['gapfilling']['gap_fill_params']['db_to_compare'] = db_to_compare
      
        if 'KEGG' in db_to_compare:
            pass
        if 'BioCyc' in db_to_compare:
            Path0 = click.prompt('Enter the path to your BioCyc TXT file containing a SmartTable with the columns \'Accession-2\' and \'Reaction of gene\'', type=click.Path(exists=True))
            Path1 = click.prompt('Enter the path to your BioCyc TXT file containing a SmartTable with all reaction relevant information', type=click.Path(exists=True))
            Path2 = click.prompt('Enter the path to your Biocyc TXT file containing a SmartTable with all metabolite relevant information', type=click.Path(exists=True))
            Path3 = click.prompt('Enter path to protein FASTA file used as input for CarveMe', type=click.Path(exists=True))
            config['gapfilling']['gap_fill_params']['biocyc_files'] = [Path0, Path1, Path2, Path3]


    # kegg pathways as groups
    kegg_pw_groups = click.prompt('Do you want to add KEGG pathways as groups to the model?', type=click.Choice(['y','n']), show_choices=True)
    config['kegg_pathway_groups'] = kegg_pw_groups

    # resolve duplicates

    reac_dups = click.prompt('Do you want to check for and/or remove duplicate reactions?', type=click.Choice(['skip','check','remove']), show_choices=True)
    config['duplicates']['reactions'] = reac_dups
    meta_dups = click.prompt('Do you want to check for and/or remove duplicate metabolites?', type=click.Choice(['skip','check','remove']), show_choices=True)
    config['duplicates']['metabolites'] = meta_dups
    unused_meta = click.prompt('Do you want to remove unused metabolites?', type=click.Choice(['y','n']), show_choices=True)
    match unused_meta:
        case 'y':
            config['duplicates']['remove_unused_metabs'] = True
        case 'n':
            config['duplicates']['remove_unused_metabs'] = False

    # BOF
    do_bofdat = click.prompt('Do you want do run BOFdat?', type=click.Choice(['y','n']), show_choices=True)
    match do_bofdat:
        case 'y':
            config['BOF']['run_bofdat'] = True
            full_genome_path = click.prompt('Please enter the path to the full genome sequence', type=click.Path(exists=True))
            config['BOF']['full_genome_sequence'] = full_genome_path
            dna_wf = click.prompt('Enter the DNA weight fraction of your organism', type=float)
            config['BOF']['dna_weight_fraction'] = dna_wf
            wf = click.prompt('Enter the wight fraction of your organsim (enzyme/ion)', type=float)
            config['BOF']['weight_fraction'] = wf
        case 'n':
            config['BOF']['run_bofdat'] = False


    # save config
    if configpath:
        pass
    else:
        configpath = Path(config['general']['dir'],'config.yaml')

    with open(configpath,'w') as outf:
        yaml.dump(config,outf,default_flow_style=False)

    return config
    
            
    

    



    

    