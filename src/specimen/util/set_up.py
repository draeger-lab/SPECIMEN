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
import warnings
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
CMPB_CONFIG_PATHS_REQUIRED = ['mediapath'] #: :meta:
CMPB_CONFIG_PATHS_OPTIONAL = ['modelpath','full_genome_sequence','gff', 'protein_fasta',
                              'gene-table','reacs-table','gff','swissprot-dmnd',
                              'swissprot-mapping'] # :meta:
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
# @DEPRECATE: MNX now coverend in refinegems - change extension after cleaning gapfill module in refinegems
def download_mnx(dir:str='MetaNetX', chunk_size:int=1024):
    """Download the data needed from the MetaNetX database.
    
    .. warning::
        This function will be deprecated soon.

    Args:
        - dir (str, optional): 
            Directory to write the downloaded files to. 
            Defaults to 'MetaNetX'.
        - chunk_size (int, optional): 
            Size of the chunk of data that is loaded into memory during download. 
            Defaults to 1024.
    """
    mes = f'This function will be deprecated soon. Functionality has been moved to refineGEMs.'
    warnings.warn(mes, FutureWarning)

    for mnx_name,mnx_url in MNX_URL_DICT.items():
        r = requests.get(mnx_url, stream=True)
        total = int(r.headers.get('content-length', 0))
        with open(Path(dir,mnx_name), 'wb') as f, tqdm(desc=str(Path(dir,mnx_name)), total=total,
            unit='iB', unit_scale=True, unit_divisor=1024,) as bar:
            for data in r.iter_content(chunk_size=chunk_size):
                size = f.write(data)
                bar.update(size)


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
            #    add part for the cmpb pipeline
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

# @TODO improve & cleanup
def validate_config(userc:str, pipeline:Literal['hqtb','cmpb']='hqtb') -> dict:
    """Validate a user hqtb config file for use in the pipeline.

    .. note::
    
        Currently not everything is checked, mainly the needed files are.

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
    
    def dict_recursive_overwrite(dictA:dict, key:str=None) -> dict:
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
                raise TypeError(F'Missing a required argument in the config file ({key}).')
            elif dictA == 'USER':
                # @TODO 
                mes = F'Keyword USER detected in config ({key}). Either due to skipped options or missing required information.\nReminder: This may lead to downstream problems.'
                logging.warning(mes)
                return None
            else:
                return dictA

        for key in dictA.keys():
            dictA[key] = dict_recursive_overwrite(dictA[key],key)
        return dictA


    # @TODO: extent - include more checks
    def dict_recursive_check(dictA:dict, key:str=None, 
                             pipeline:Literal['hqtb','cmpb']='hqtb'):
        """Helper-function for :py:func:`~specimen.util.set_up.validate_config` 
        to check if a configuration is valid to run the high-quality template based pipeline.

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
                if isinstance(dictA,list):
                    for entry in dictA:
                        if os.path.isfile(entry):
                            continue
                        else:
                            raise FileNotFoundError(F'Path does not exist: {dictA}')
                elif dictA and os.path.isfile(dictA):
                    return
                else: 
                    raise FileNotFoundError(F'Path does not exist: {dictA}')
            # optional file paths
            elif key in PIPELINE_PATHS_OPTIONAL[pipeline]:
                if isinstance(dictA,str):
                    if os.path.isfile(dictA):
                        return
                    elif not os.path.isfile(dictA): 
                        mes = F'Path does not exist: {dictA}. \nReminder: It is optional, but it may lead to downstream problems.'
                        logging.warning(mes)
                        pass
                    else:
                        raise FileNotFoundError(F'Path does not exist: {dictA}')
                if isinstance(dictA,list):
                    for entry in dictA:
                        if entry and os.path.isfile(entry):
                            return
                        elif not os.path.isfile(entry):
                            mes = F'Path does not exist: {entry}. \nReminder: It is optional, but it may lead to downstream problems.'
                            logging.warning(mes)
                            pass
                        else: 
                            raise FileNotFoundError(F'Path does not exist: {entry}')
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

    if combined_config['general']['modelname'] is None and (combined_config['general']['authorinitials'] is None or combined_config['general']['organism'] is None or combined_config['general']['strainid'] is None):
        raise ValueError(f'Either the model name or all of the following parameters must be stated: authorinitials, organism and strainID')

    return combined_config


# cmpb
# ----

def save_cmpb_user_input(configpath:Union[str,None]=None) -> dict:
    """Guide the user step by step through the creation of the configuration for a cmpb pipeline run
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
    config['input']['mediapath'] = click.prompt('Enter the path to a media configuration file for growth simulation', type=click.Path(exists=True))

    # general options
    # ---------------
    print('------------')
    print('General options')
    print('------------')

    # output directory
    config['general']['dir'] = click.prompt('Enter your desired output directory path', type=click.Path())

    # name for the model
    modelname = click.prompt('Do you have a specific name for your model?', type=click.Choice(['y','n']), show_choices=True)
    match modelname:
        case 'y':
            config['carveme']['modelname'] = click.prompt('Please enter your desired name for the model', type=str)
        case 'n':
            config['general']['authorinitials'] = click.prompt('An automated name based on the pattern iOrganismStrainAuthorYear will be created. \n Please enter your intials.', type=str)
            config['general']['organism'] = click.prompt('Please enter an abbreviation for your organism.', type=str)
            config['general']['strainid'] = click.prompt('Please enter the ID for your strain.', type=str)
    
    # colour 
    set_col = click.prompt('Do you want to use the default colour map YlGn for the visualisation?', type=click.Choice(['y','n']), show_choices=True)
    match set_col:
        case 'n':
            colours = click.prompt('Enter your chosen colour scheme', type=str)
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
    run_memote = click.prompt('Do you want to run memote after each step?', type=click.Choice(['y','n']), show_choices=True)
    match run_memote:
        case 'y':
            config['general']['memote_always_on'] = True
        case 'n':
            config['general']['memote_always_on'] = False
    
    # run stats always y/n
    models_stats = click.prompt('Do you want to run stats after each step?', type=click.Choice(['y','n']), show_choices=True)
    match models_stats:
        case 'y':
            config['general']['stats_always_on'] = True
        case 'n':
            config['general']['stats_always_on'] = False

    # some additional, sometimes required, sometimes optional files
    refseq = click.prompt('If you want to run a gap analysis with KEGG or have a CarveMe model, please enter the path to your refseq gff file', type=click.Path())
    config['general']['gff'] = refseq

    kegg_org_id = click.prompt('If you want to run a gap analysis with KEGG, please enter the KEGG organism ID')
    config['general']['kegg_organism_id'] = kegg_org_id
    
    protein_fasta = click.prompt('If you want to use CarveMe or GeneGapFiller, please enter the path to your protein fasta file', type=click.Path())
    config['general']['protein_fasta'] = protein_fasta

    # tech resources
    # --------------
    email = click.prompt('Enter the e-mail that will be used for Entrez')
    config['tech-resources']['email'] = email

    set_threads = click.prompt('The default number of threads available for tools like DIAMOND is 2. Do you want to change that?', type=click.Choice(['y','n']), show_choices=True)
    match set_threads:
        case 'y':
            threads = click.prompt('Enter the number of threads available for tools like DIAMOND', type=int)
            config['tech-resources']['threads'] = threads

    # part-specific 
    # -------------
    print('------------')
    print('Part-specific options')
    print('------------')

    # CarveMe
    carve = click.prompt('Do you want to build a model using CarveMe?', type=click.Choice(['y','n']), show_choices=True)
    match carve:
        case 'y':
            if config['general']['protein_fasta'] is None:
                protein_fasta = click.prompt('Enter the path to your protein fasta file', type=click.Path(exists=True))
                config['general']['protein_fasta'] = protein_fasta
            gram = click.prompt('Do you want to use a template specialized for gram-positive or gram-negative bacteria?', type=click.Choice(['grampos','gramneg','None']), show_choices=True)
            config['carveme']['gram'] = gram
        case 'n':
            if config['input']['modelpath'] is None:
                model = click.prompt('Please choose between an existing model or building a model with CarveMe. To run the CMPB workflow, you need a model.', type=click.Choice(['modelpath','CarveMe']), show_choices=True)
                match model:
                    case 'modelpath':
                        modelpath = click.prompt('Enter the path to an existing model', type=click.Path(exists=True))
                        config['input']['modelpath'] = modelpath
                    case 'CarveMe':
                        carveme = click.prompt('Enter the path to a protein fasta file', type=click.Path(exists=True))
                        config['general']['protein_fasta'] = carveme
    
    # model polish
    carveme = click.prompt('Is your draft model CarveMe-based?', type=click.Choice(['y','n']), show_choices=True)
    if carveme == 'y':
        labs = click.prompt('Do you have a strain without any database information?', type=click.Choice(['y','n']), show_choices=True)
        labs = True if labs == 'y' else False
        config['cm-polish']['is_lab_strain'] = labs
    
    # gapfilling
    gap_analysis = click.prompt('Do you want to run a gap analysis?', type=click.Choice(['y','n']), show_choices=True) 

    if gap_analysis == 'y':
        idprefix = click.prompt('Enter a prefix to be used of IDs for the namespace do not exist')
        config['gapfilling']['idprefix'] = idprefix
        formula_check = click.prompt('Enter the parameter for checking the metabolite formula before adding them to the model', type=click.Choice(['none','strict','existence','wildcard']), show_choices=True)
        config['gapfilling']['formula-check'] = formula_check
        exclude_dna = click.prompt('Do you want to exlude reactions containing \'DNA\' in their name?', type=click.Choice(['y','n']), show_choices=True)
        config['gapfilling']['exclude-dna'] = exclude_dna
        exclude_rna = click.prompt('Do you want to exlude reactions containing \'RNA\' in their name?', type=click.Choice(['y','n']), show_choices=True)
        config['gapfilling']['exclude-rna'] = exclude_rna
        config['gapfilling']['threshold_add_reacs'] = click.prompt('Enter the threshold for adding reactions (max. allowed matches of an EC number).', type=int, default=5)
        
        algorithm = click.prompt('Which algorithm do you want to use for gapfilling?', type=click.Choice(['KEGGapFiller','BioCycGapFiller','GeneGapFiller']), show_choices=True)
        another_gapfiller = True
        while another_gapfiller:
            match algorithm:
                case 'KEGGapFiller':
                    config['gapfilling']['KEGGapFiller'] = True

                    if config['general']['kegg_organism_id'] is None:
                        kegg_org_id = click.prompt('Enter the KEGG organism id')
                        config['general']['kegg_organism_id'] = kegg_org_id
                case 'BioCycGapFiller':
                    config['gapfilling']['BioCycGapFiller'] = True

                    gene_table = click.prompt('Enter the path to a gene smart table from BioCyc', type=click.Path(exists=True))
                    config['gapfilling']['BioCycGapFiller parameters']['gene-table'] = gene_table
                    reacs_table = click.prompt('Enter the path to a reactions smart table from BioCyc', type=click.Path(exists=True))
                    config['gapfilling']['BioCycGapFiller parameters']['reacs-table'] = reacs_table
                    gff = click.prompt('Enter the path to a GFF file of the genome of the model', type=click.Path(exists=True))
                    config['gapfilling']['BioCycGapFiller parameters']['gff'] = gff
                case 'GeneGapFiller':
                    config['gapfilling']['GeneGapFiller'] = True

                    gff = click.prompt('Enter the path to a GFF file of the genome of the model', type=click.Path(exists=True))
                    config['gapfilling']['GeneGapFiller parameters']['gff'] = gff
                    swissprot_dmnd = click.prompt('Enter the path to the SwissProt DIAMOND database file', type=click.Path(exists=True))
                    config['gapfilling']['GeneGapFiller parameters']['swissprot-dmnd'] = swissprot_dmnd
                    swissprot_mapping = click.prompt('Enter the path to the SwissProt mapping file', type=click.Path(exists=True))
                    config['gapfilling']['GeneGapFiller parameters']['swissprot-mapping'] = swissprot_mapping
                    check_NCBI = click.prompt('Do you want to enable checking NCBI accession numbers for EC numbers?', type=click.Choice(['y','n']), show_choices=True)
                    check_NCBI = True if check_NCBI == 'y' else False
                    config['gapfilling']['GeneGapFiller parameters']['check-NCBI'] = check_NCBI
                    sensitivity = click.prompt('Enter the sensitivity option for the DIAMOND run', type=click.Choice(['fast','mid-sensitive','sensitive','more-sensitive','very-sensitive','ultra-sensitive']), show_choices=True)
                    config['gapfilling']['GeneGapFiller parameters']['sensitivity'] = sensitivity
                    coverage = click.prompt('Enter the coverage for DIAMOND', type=float)
                    config['gapfilling']['GeneGapFiller parameters']['coverage'] = coverage
                    percentage_identity = click.prompt('Enter the percentage identity threshold value for accepting matches', type=float)
                    config['gapfilling']['GeneGapFiller parameters']['percentage identity'] = percentage_identity
            another_gapfiller = click.prompt('Do you want to use another algorithm for gapfilling?', type=click.Choice(['y','n']), show_choices=True)
            another_gapfiller = True if another_gapfiller == 'y' else False
            if another_gapfiller:
                algorithm = click.prompt('Which algorithm do you want to use for gapfilling?', type=click.Choice(['KEGGapFiller','BioCycGapFiller','GeneGapFiller']), show_choices=True)

    # ModelPolisher
    modelpolisher = click.prompt('Do you want to run ModelPolisher?', type=click.Choice(['y','n']), show_choices=True)
    match modelpolisher:
        case 'y':
            config['modelpolisher'] = True
            allow_model_to_be_saved_on_server = click.prompt('Do you want to allow the model to be saved on the server?', type=click.Choice(['y','n']), show_choices=True)
            allow_model_to_be_saved_on_server = True if allow_model_to_be_saved_on_server == 'y' else False
            config['mp']['allow-model-to-be-saved-on-server'] = allow_model_to_be_saved_on_server

            dont_fix = click.prompt('Do you want to fix the model? Unset default values will be set, if they are mandatory.', type=click.Choice(['y','n']), show_choices=True)
            dont_fix = False if dont_fix == 'y' else True
            config['mp']['fixing']['dont-fix'] = dont_fix
            
            annotate_with_bigg = click.prompt('Do you want to annotate with BiGG?', type=click.Choice(['y','n']), show_choices=True)
            annotate_with_bigg = True if annotate_with_bigg == 'y' else False
            config['mp']['annotation']['bigg']['annotate-with-bigg'] = annotate_with_bigg
            include_any_uri = click.prompt('Do you want to include annotation that are not MIRIAM-compliant?', type=click.Choice(['y','n']), show_choices=True)
            include_any_uri = True if include_any_uri == 'y' else False
            config['mp']['annotation']['bigg']['include-any-uri'] = include_any_uri
        case 'n':
            config['modelpolisher'] = False

    # kegg pathways as groups
    kegg_pw_groups = click.prompt('Do you want to add KEGG pathways as groups to the model?', type=click.Choice(['y','n']), show_choices=True)
    match kegg_pw_groups:
        case 'y':
            config['kegg_pathway_groups'] = True
        case 'n':
            config['kegg_pathway_groups'] = False

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
            
    # handling EGCs
    egc_solver = click.prompt('Choose a solver (or none) for handling energy generating cycles.', type=click.Choice(['none','greedy']), show_choices=True)
    if egc_solver == 'none':
        egc_solver = None

    # BOF
    do_bofdat = click.prompt('Do you want do run BOFdat?', type=click.Choice(['y','n']), show_choices=True)
    match do_bofdat:
        case 'y':
            config['BOF']['run_bofdat'] = True
            full_genome_path = click.prompt('Please enter the path to the full genome sequence', type=click.Path(exists=True))
            config['BOF']['full_genome_sequence'] = full_genome_path
            dna_wf = click.prompt('Enter the DNA weight fraction of your organism', type=float)
            config['BOF']['dna_weight_fraction'] = dna_wf
            wf = click.prompt('Enter the weight fraction of your organism (enzyme/ion)', type=float)
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
    
            
    

    



    

    