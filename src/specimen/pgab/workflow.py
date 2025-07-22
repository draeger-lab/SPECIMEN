#!/usr/bin/env python3
"""Functions to run the workflow based on the Prokaryotic Genome Annotation Pipeline (PGAP). """

from typing import Union, Literal
from pathlib import Path
import yaml
import logging
import os
import subprocess
import xml.etree.ElementTree as ET
import pandas as pd
from Bio import SeqIO
from tqdm.auto import tqdm
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import datetime
from datetime import date
from importlib.resources import files
from ..cmpb.workflow import run as run_cmpb
from ..hqtb.workflow import run as run_hqtb
import warnings
import requests
import gzip
import shutil

# config keys for pipeline files
HQTB_CONFIG_PATH_OPTIONAL = [
    "media_gap",
    "ncbi_map",
    "biocyc",
    "universal",
    "pan-core",
    "fasta",
    "gff",
    "dmnd-database",
    "database-mapping",
]  #: :meta:
HQTB_CONFIG_PATH_REQUIRED = [
    "annotated_genome",
    "full_sequence",
    "model",
    "diamond",
    "media_analysis",
]  #: :meta:
CMPB_CONFIG_PATHS_REQUIRED = ["mediapath"]  #: :meta:
CMPB_CONFIG_PATHS_OPTIONAL = [
    "modelpath",
    "full_genome_sequence",
    "gff",
    "protein_fasta",
    "gene-table",
    "reacs-table",
    "gff",
    "dmnd-database",
    "database-mapping",
    "reaction_direction",
]  # :meta:
PIPELINE_PATHS_OPTIONAL = {
    "hqtb": HQTB_CONFIG_PATH_OPTIONAL,
    "cmpb": CMPB_CONFIG_PATHS_OPTIONAL,
    "pgab": ["configpath","db"]
}  #: :meta:
PIPELINE_PATHS_REQUIRED = {
    "hqtb": HQTB_CONFIG_PATH_REQUIRED,
    "cmpb": CMPB_CONFIG_PATHS_REQUIRED,
    "pgab": ["fasta"]
}  #: :meta:
# config keys for pipelines directories
PIPELINE_DIR_PATHS = ["dir"]

# TODO: look into warnings instead of print/return
# TODO: improve variable names
# TODO: outputs for the current step
# TODO: log-file?
# TODO: information for each of the methods
def run(configpath: Union[str, None] = None):
    # load config
    # -----------
    if not configpath:
        # TODO: user input method for PGAB
        config = save_pgab_user_input()
    else:
        # TODO: adapt in case of changes, works for now
        config = validate_config(configpath, "pgab")

    # create directory structure
    # --------------------------
    # TODO: add required folders
    # build_data_directories("pgab", config['general']['dir'])

    dir = str(Path(config['general']['dir'], 'pgab_out'))

    # set up model name
    # -----------------
    if (
        config["general"]["authorinitials"] is not None
        and config["general"]["organism"] is not None
        and config["general"]["strainid"] is not None
    ):
        modelname = (
            "i"
            + config["general"]["organism"]
            + str(config["general"]["strainid"])
            + config["general"]["authorinitials"]
            + str(date.today().year).removeprefix("20")
        )
    elif config["general"]["modelname"] is not None:
        modelname = config["general"]["modelname"]
    else:
        print(
            "No values given for the standard name for a model. Default name will be used."
        )
        modelname = "model_" + str(date.today().year).removeprefix("20")

    # run PGAP
    # --------
    # run_pgap(config['pgap'], str(Path(dir,'pgap')))

    match_100 = None
    match config['pgap']['tax-check']:
        case 'only':
            print('The taxonomy check with PGAP is finished. The results can be viewed in XML- or txt-format.')
            result = parse_tax_check(str(Path(dir,'pgap','output','ani-tax-report.xml')), type='only')
            print('As a further step, running the HQTB pipeline with one of these strains as a template is recommended:')
            print(result.to_string(index=False))
            return
        case 'continue':
            taxonomy, match_100 = parse_tax_check(str(Path(dir,'pgap','output','ani-tax-report.xml')), config['pgap']['tax_cutoff'])
            # transform dict into a list for DIAMOND
            if len(taxonomy) > 1:
                taxid = [taxonomy['best'], taxonomy['predicted']]
            else: taxid = [taxonomy['predicted']]
            
            # TODO: default cutoff?
            if taxid:
                if config['pgap']['skip_diamond']=='yes':
                    print(f'The taxid(s) which will be used for the taxonomy specific DIAMOND run are {taxid}.')
            else:
                print('There was no taxid found above the cutoff. You will have to decide on a taxid manually.')
                return
        case 'none':
            if config['pgap']['taxid']:
                taxid = str(config['pgap']['taxid']).split(',')
                taxid = [t.strip() for t in taxid]
            else: taxid = None
        case _:
            raise ValueError(f"Unknown input for taxonomy check: {config['pgap']['tax-check']}")

    # option to skip DIAMOND if a 100% match is found
    if match_100 and config['pgap']['skip_diamond']!='no':
        # TODO: maybe revise the structure?
        print(f'DIAMOND will be skipped because of a found organism with 100% identity. As a next step, the {config['polish']['next_step']} workflow will be started.')

        # converting the Genbank accession to a RefSeq accession
        with open('result.txt', 'w') as f: # TODO: merge both open() into one?
            result = subprocess.run(['./datasets', 'summary', 'genome', 'accession', match_100['accession']], stdout=f)
        with open('result.txt', 'r') as f:
            result = f.readline()
        if 'GCF' in result:
            result = result.split(',')
            for item in result:
                if 'GCF' in item and 'paired_accession' in item:
                    refseq = item.split('\"')[3]
            
            # run next workflow with the RefSeq accession
            # QUESTION: why is the directory for HQTB ./specimen_run/?
            match config['polish']['next_step']:
                case 'cmpb':
                    polish_config = adapt_config(config['polish'], dir, modelname, refseq=refseq)
                    run_cmpb(polish_config)
                    # print('work in progress...')
                case 'hqtb':
                    accession = refseq.split('_')[1].split('.')[0]
                    numbers = [accession[i:i+3] for i in range(0, len(accession), 3)]
                    accession = f"{refseq}_ASM{accession.lstrip('0').rstrip('5')}v{refseq.split('.')[1]}"
                    response = requests.get(f"https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/{numbers[0]}/{numbers[1]}/{numbers[2]}/{accession}/{accession}_translated_cds.faa.gz")
                    if response.status_code == 200:
                        with open(str(Path(dir,"annot_translated_cds.faa.gz")), mode="wb") as file:
                            file.write(response.content)
                        with gzip.open(str(Path(dir,"annot_translated_cds.faa.gz")), 'rb') as file_in:
                            with open(str(Path(dir,"annot_translated_cds.faa")), 'wb') as file_out:
                                shutil.copyfileobj(file_in, file_out)
                        config['polish']['fasta'] = str(Path(dir,"annot_translated_cds.faa"))
                    else:
                        print(f"There was no protein FASTA file found for {refseq} the automated way. To run HQTB, you will have to specify the files yourself.")
                        return
                    response = requests.get(f"https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/{numbers[0]}/{numbers[1]}/{numbers[2]}/{accession}/{accession}_genomic.gff.gz")
                    if response.status_code == 200:
                        with open(str(Path(dir,"annot.gff.gz")), mode="wb") as file:
                            file.write(response.content)
                        with gzip.open(str(Path(dir,"annot.gff.gz")), 'rb') as file_in:
                            with open(str(Path(dir,"annot.gff")), 'wb') as file_out:
                                shutil.copyfileobj(file_in, file_out)
                        config['polish']['gff'] = str(Path(dir,"annot.gff"))
                    else:
                        print(f"There was no GFF file found for {refseq} the automated way. To run HQTB, you will have to specify the files yourself.")
                        return
                    polish_config = adapt_config(config['polish'], dir, modelname, refseq, genome=config['pgap']['generic']['fasta']['location'])
                    run_hqtb(polish_config)
                case _:
                    raise ValueError(f"Unknown input for next polishing step: {config['polish']['next_step']}")
            return
        else: 
            print('There was no matching RefSeq-ID found and CarveMe could not be performed.')
            if config['pgap']['skip_diamond']=='maybe':
                print(f'Because of this, DIAMOND will run with the following taxonomy id(s): {taxid}')
            else: return

    # first save of FASTA with "right" file name
    fasta_nr = str(Path(dir,'pgap','output','annot_translated_cds.faa'))
    fasta_sp = str(Path(dir,'pgap','output','annot_translated_cds.faa'))
    mapping_nr = []        # mapping table between locus tag and protein id (just one, gets appended) [1. output]
    duplicates_nr = {}     # removed duplicates (list for each run) [3. output]
    below_cutoff_nr = {}   # bad identity (list for each run) [4. output]
    statistics_nr = {}     # statistics for every DIAMOND run (dictionary) [2. output]
    mapping_sp = []        # mapping table between locus tag and protein id (just one, gets appended) [1. output]
    duplicates_sp = {}     # removed duplicates (list for each run) [3. output]
    below_cutoff_sp = {}   # bad identity (list for each run) [4. output]
    statistics_sp = {}     # statistics for every DIAMOND run (dictionary) [2. output]

    if taxid is not None:
        for id in taxid:
            run_id = []
            run_id.append(id)
            
            # nr database
            if config['nr_blast']['run']:
                print(datetime.datetime.now())
                print(f'DIAMOND run with taxid {id} and nr database')
                run_DIAMOND(fasta_nr, config['nr_blast']['db'], 
                            config['nr_blast']['sensitivity'], config['nr_blast']['coverage'], config['nr_blast']['threads'], 
                            str(Path(dir,'nr_blast')), run_id, f'nr_{id}')
                stats_tax, duplicates_tax, cutoff_tax = get_mapping_table(mapping_nr, str(Path(dir,'nr_blast',f'nr_{id}_{config['nr_blast']['sensitivity']}.tsv')), fasta_nr, dir, config['nr_blast']['pid'], 'nr')
                statistics_nr[id] = stats_tax
                duplicates_nr[id] = duplicates_tax
                below_cutoff_nr[id] = cutoff_tax

            fasta_nr = str(Path(dir,'DIAMOND_input_nr.faa'))
            
            # SwissProt
            if config['swissprot_blast']['run']:
                print(datetime.datetime.now())
                print(f'DIAMOND run with taxid {id} and SwissProt database')
                run_DIAMOND(fasta_sp, config['swissprot_blast']['db'], 
                            config['swissprot_blast']['sensitivity'], config['swissprot_blast']['coverage'], config['swissprot_blast']['threads'], 
                            str(Path(dir,'swissprot_blast')), run_id, f'swissprot_{id}')
                stats_tax, duplicates_tax, cutoff_tax = get_mapping_table(mapping_sp, str(Path(dir,'swissprot_blast',f'swissprot_{id}_{config['swissprot_blast']['sensitivity']}.tsv')), fasta_sp, dir, config['swissprot_blast']['pid'], 'swissprot')
                statistics_sp[id] = stats_tax
                duplicates_sp[id] = duplicates_tax
                below_cutoff_sp[id] = cutoff_tax
            
            fasta_sp = str(Path(dir,'DIAMOND_input_swissprot.faa'))

    # run DIAMOND without taxonomy sensitivity
    if config['nr_blast']['run']:
        print(datetime.datetime.now())
        print('DIAMOND run without taxonomy restriction and nr database')
        run_DIAMOND(fasta_nr, config['nr_blast']['db'], 
                    config['nr_blast']['sensitivity'], config['nr_blast']['coverage'], config['nr_blast']['threads'], 
                    str(Path(dir,'nr_blast')), outname='nr_not_taxrestricted')
        print(datetime.datetime.now())

        # Expand mapping table
        stats_notax, duplicates_notax, cutoff_notax = get_mapping_table(mapping_nr, str(Path(dir,'nr_blast','nr_not_taxrestricted_'+config['nr_blast']['sensitivity']+'.tsv')), fasta_nr, dir, config['nr_blast']['pid'], 'nr')
        statistics_nr['no_tax_restriction'] = stats_notax
        duplicates_nr['no_tax_restriction'] = duplicates_notax
        below_cutoff_nr['no_tax_restriction'] = cutoff_notax
        
        # DIAMOND report
        report = DiamondReport(statistics_nr, duplicates_nr, below_cutoff_nr)
        report.save(str(Path(dir,'nr_blast')))
        
    if config['swissprot_blast']['run']:
        print(datetime.datetime.now())
        print('DIAMOND run without taxonomy restriction and SwissProt database')
        run_DIAMOND(fasta_sp, config['swissprot_blast']['db'], 
                    config['swissprot_blast']['sensitivity'], config['swissprot_blast']['coverage'], config['swissprot_blast']['threads'], 
                    str(Path(dir,'swissprot_blast')), outname='swissprot_not_taxrestricted')
        print(datetime.datetime.now())

        # Expand mapping table
        stats_notax, duplicates_notax, cutoff_notax = get_mapping_table(mapping_sp, str(Path(dir,'swissprot_blast','swissprot_not_taxrestricted_'+config['swissprot_blast']['sensitivity']+'.tsv')), fasta_sp, dir, config['swissprot_blast']['pid'], 'swissprot')
        statistics_sp['no_tax_restriction'] = stats_notax
        duplicates_sp['no_tax_restriction'] = duplicates_notax
        below_cutoff_sp['no_tax_restriction'] = cutoff_notax
        
        # DIAMOND report
        report = DiamondReport(statistics_sp, duplicates_sp, below_cutoff_sp)
        report.save(str(Path(dir,'swissprot_blast')))
        
    # save mapping tables (even if empty)
    mapping_nr = pd.DataFrame(mapping_nr, columns=['model_id','NCBI'])
    mapping_nr.to_csv(str(Path(dir, 'mapping_table_nr.csv')), header=True, index=False, sep='\t')
    # QUESTION: UniProt instead of DECLASSIFIED?
    mapping_sp = pd.DataFrame(mapping_sp, columns=['model_id','UNIPROT'])
    mapping_sp.to_csv(str(Path(dir, 'mapping_table_swissprot.csv')), header=True, index=False, sep='\t')

    # TODO: extra BLAST manually, user input?

    # QUESTION: Why doesnt CarveMe add specific proteins into the model?
        
    polish_config = adapt_config(config['polish'], dir, modelname, genome=config['pgap']['generic']['fasta']['location'])
    match config['polish']['next_step']:
        case 'cmpb':
            run_cmpb(polish_config)
            # print('work in progress...')
        case 'hqtb':
            # TODO: automated way to find a template model?
            run_hqtb(polish_config)
            # print('work in progress...')
        case _:
            raise ValueError(f"Unknown input for next polishing step: {config['polish']['next_step']}")

    # QUESTION: What exactly belongs in the header of the mapping table? 'NCBI' for protein_tag and 'UNIPROT' for UniProt?
    
def save_pgab_user_input() -> dict[str, str]:
    """method adapted from SPECIMEN.util.set_up
    !under construction!
    """

def validate_config(userc: str, pipeline: Literal["hqtb", "cmpb", "pgab"] = "hqtb") -> dict:
    """method adapted from SPECIMEN.util.set_up"""

    def dict_recursive_combine(dictA: dict, dictB: dict) -> dict:
        if not isinstance(dictB, dict):
            return dictB
        for key in dictA.keys():
            if key in dictB.keys():
                dictA[key] = dict_recursive_combine(dictA[key], dictB[key])
        return dictA

    def dict_recursive_overwrite(dictA: dict, key: str = None) -> dict:
        if not isinstance(dictA, dict):
            # check for missing input
            if dictA == "__USER__":
                raise TypeError(
                    f"Missing a required argument in the config file ({key})."
                )
            elif dictA == "USER":
                mes = f"Keyword USER detected in config ({key}). Either due to skipped options or missing required information.\nReminder: This may lead to downstream problems."
                logging.warning(mes)
                return None
            else:
                return dictA

        for key in dictA.keys():
            dictA[key] = dict_recursive_overwrite(dictA[key], key)
        return dictA

    def dict_recursive_check(dictA: dict, key: str = None, pipeline: Literal["hqtb", "cmpb", "pgab"] = "hqtb"):
        if not isinstance(dictA, dict):
            # required file paths
            if key in PIPELINE_PATHS_REQUIRED[pipeline]:
                if isinstance(dictA, list):
                    for entry in dictA:
                        if os.path.isfile(entry):
                            continue
                        else:
                            raise FileNotFoundError(f"Path does not exist: {dictA}")
                elif dictA and os.path.isfile(dictA):
                    return
                else:
                    raise FileNotFoundError(f"Path does not exist: {dictA}")
            # optional file paths
            elif key in PIPELINE_PATHS_OPTIONAL[pipeline]:
                if isinstance(dictA, str):
                    if os.path.isfile(dictA):
                        return
                    elif not os.path.isfile(dictA):
                        mes = f"Path does not exist: {dictA}. \nReminder: It is optional, but it may lead to downstream problems."
                        logging.warning(mes)
                        pass
                    else:
                        raise FileNotFoundError(f"Path does not exist: {dictA}")
                if isinstance(dictA, list):
                    for entry in dictA:
                        if entry and os.path.isfile(entry):
                            return
                        elif not os.path.isfile(entry):
                            mes = f"Path does not exist: {entry}. \nReminder: It is optional, but it may lead to downstream problems."
                            logging.warning(mes)
                            pass
                        else:
                            raise FileNotFoundError(f"Path does not exist: {entry}")
            elif key in PIPELINE_DIR_PATHS:
                if dictA and os.path.exists(dictA):
                    return
                else:
                    raise FileNotFoundError(f"Directory does not exist: {dictA}")
            # not found or missing
            else:
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
        case "hqtb":
            defaultc_path = files("specimen.data.config").joinpath("hqtb_config_default.yaml")
        case "cmpb":
            defaultc_path = files("specimen.data.config").joinpath("cmpb_config.yaml")
        case "pgab":
            defaultc_path = files("specimen.data.config").joinpath("pgab_config.yaml")
            # defaultc_path = "/vol/pgap_data/pgab_config.yaml"
        case _:
            raise ValueError(f"Unknown input for pipeline: {pipeline}")

    with open(defaultc_path, "r") as cfg_def, open(userc, "r") as cfg_usr:
        config_d = yaml.load(cfg_def, Loader=yaml.loader.FullLoader)
        config_u = yaml.load(cfg_usr, Loader=yaml.loader.FullLoader)

    # combine
    combined_config = dict_recursive_combine(config_d, config_u)

    # overwrite __USER__ and USER
    combined_config = dict_recursive_overwrite(combined_config)

    # check for missing or problematic values
    # special case for HQTB pipeline with relative paths
    if (
        "data" in combined_config.keys()
        and "data_direc" in combined_config["data"].keys()
        and combined_config["data"]["data_direc"]
    ):
        if os.path.isdir(combined_config["data"]["data_direc"]):
            for key in combined_config["data"]:
                if combined_config["data"][key] and key != "data_direc":
                    combined_config["data"][key] = (
                        combined_config["data"]["data_direc"]
                        + combined_config["data"][key]
                    )
            dict_recursive_check(combined_config, key=None, pipeline=pipeline)
        else:
            raise FileNotFoundError(
                "Directory set for config:data:data_direc does not exist."
            )
    # normal recursion for validation
    else:
        dict_recursive_check(combined_config, key=None, pipeline=pipeline)

    if combined_config["general"]["modelname"] is None and (
        combined_config["general"]["authorinitials"] is None
        or combined_config["general"]["organism"] is None
        or combined_config["general"]["strainid"] is None
    ):
        raise ValueError(
            f"Either the model name or all of the following parameters must be stated: authorinitials, organism and strainID"
        )

    return combined_config

def build_data_directories(pipeline: Literal["hqtb", "high-quality template based", "cmpb", "carveme modelpolisher based", "pgab", "PGAP based"], dir: str):
    """method adapted from SPECIMEN.util.set_up"""

    match pipeline:
        # HQTB setup
        case "hqtb" | "high-quality template based":
            # create the data directory structure
            print("Creating directory structure...")
            DATA_DIRECTORIES = [
                "annotated_genomes",
                "BioCyc",
                "RefSeqs",
                "medium",
                "pan-core-models",
                "template-models",
                "universal-models",
            ]
            for sub_dir in DATA_DIRECTORIES:
                new_dir = Path(dir, sub_dir)
                try:
                    Path(new_dir).mkdir(parents=True, exist_ok=False)
                    print(f"Creating new directory {new_dir}")
                except FileExistsError:
                    print(f"Directory {new_dir} already exists.")

        # CMPB output
        case "cmpb" | "carveme modelpolisher based":
            Path(dir, "cmpb_out").mkdir(parents=True, exist_ok=False)  # cmpb_out
            Path(dir, "cmpb_out", "models").mkdir(
                parents=True, exist_ok=False
            )  #   |- models
            Path(dir, "cmpb_out", "logs").mkdir(
                parents=True, exist_ok=False
            )  #   |- logs
            Path(dir, "cmpb_out", "misc").mkdir(
                parents=True, exist_ok=False
            )  #   |- misc
            Path(dir, "cmpb_out", "misc", "memote").mkdir(
                parents=True, exist_ok=False
            )  #      |- memote
            Path(dir, "cmpb_out", "misc", "mcc").mkdir(
                parents=True, exist_ok=False
            )  #      |- mcc
            Path(dir, "cmpb_out", "misc", "gapfill").mkdir(
                parents=True, exist_ok=False
            )  #      |- gapfill
            Path(dir, "cmpb_out", "misc", "growth").mkdir(
                parents=True, exist_ok=False
            )  #      |- growth
            Path(dir, "cmpb_out", "misc", "stats").mkdir(
                parents=True, exist_ok=False
            )  #      |- stats
            Path(dir, "cmpb_out", "misc", "modelpolisher").mkdir(
                parents=True, exist_ok=False
            )  #      |- modelpolisher
            Path(dir, "cmpb_out", "misc", "kegg_pathway").mkdir(
                parents=True, exist_ok=False
            )  #      |- kegg_pathways
            Path(dir, "cmpb_out", "misc", "auxotrophy").mkdir(
                parents=True, exist_ok=False
            )  #      |- auxothrophy

        case "pgab" | "PGAP based":
            print("Creating directory structure...")
            Path(dir, "pgab_out").mkdir(parents=True, exist_ok=False)  # pgab_out
            Path(dir, "pgab_out", "pgap").mkdir(
                parents=True, exist_ok=False
            )  #   |- yaml files + PGAP output
            Path(dir, "pgab_out", "nr_blast").mkdir(
                parents=True, exist_ok=False
            )  #   |- DIAMOND nr output
            Path(dir, "pgab_out", "swissprot_blast").mkdir(
                parents=True, exist_ok=False
            )  #   |- DIAMOND swissprot output

        # default case
        case _:
            message = f"Unknown input for parameter pipeline: {pipeline}"
            raise ValueError(message)

def run_pgap(config: dict, dir: str):
    logfile = str(Path(dir,'log_pgap.txt'))

    prepare_pgap_input(config['generic'], config['metadata'], dir)

    match config['tax-check']:
        # BUG: does not work without --ignore-all-errors
        case 'only': # user searches an organism, no PGAP
            completed_pgap = subprocess.run(['./pgap.py', '-r', '--taxcheck-only', '--no-self-update', '--ignore-all-errors', '-o', str(Path(dir, 'output')), str(Path(dir, 'input.yaml'))],
                shell=False, stderr=subprocess.PIPE, text=True)
        case 'continue': # user does taxcheck and PGAP
            completed_pgap = subprocess.run(['./pgap.py', '-r', '--taxcheck', '--no-self-update', '--ignore-all-errors', '-o', str(Path(dir, 'output')), str(Path(dir, 'input.yaml'))],
                shell=False, stderr=subprocess.PIPE, text=True)
        case 'none': # user already has an organism
            completed_pgap = subprocess.run(['./pgap.py', '-r', '--no-self-update', '--ignore-all-errors', '-o', str(Path(dir, 'output')), str(Path(dir, 'input.yaml'))],
                shell=False, stderr=subprocess.PIPE, text=True)
        case _:
            raise ValueError(f"Unknown input for taxonomy check: {config['pgap']['tax-check']}")
    with open(logfile, "a") as f:
        f.write(completed_pgap.stderr)

def prepare_pgap_input(generic: dict, metadata: str, dir: str):
    subprocess.run(['cp', generic['fasta']['location'], dir])

    with open(Path(dir,'submol.yaml'), 'w') as f:
        yaml.dump(metadata, f)

    generic.update({'submol': {'class': 'File', 'location': str(Path(dir,'submol.yaml'))}})
    generic['fasta']['location'] = generic['fasta']['location'].split('/')[-1]

    with open(Path(dir,'input.yaml'), 'w') as f:
        yaml.dump(generic, f)

def parse_tax_check(file: Path, cutoff: float=90.0, type: Literal['only','continue'] = 'continue') -> tuple[dict[str, str], dict[str, str]|None]:
    taxonomy = {}
    results = []
    child_100 = {}

    data = ET.parse(file)
    root = data.getroot()

    for child in root[1]:
        if child.tag == 'subject' :
            # Search for a 100% match
            if float(child.get('ANI')) == 100:
                print(f"An organism with a 100% identity was found. This organism is {child.get('org-name')}.")
                child_100['taxid'] = child.get('taxid')
                child_100['accession'] = child.get('asm_accession')
            # Search for matches above the cutoff
            if float(child.get('ANI')) >= cutoff:
                results.append({'organism':child.get('org-name'), 'identity':float(child.get('ANI')), 'taxid':child.get('taxid'), 'accession':child.get('asm_accession')})
    
    df = pd.DataFrame(results).sort_values(by='identity', ascending=False)
    if type == 'only':
        return df
    taxonomy['best'] = df.loc[1, 'taxid']

    # Search for the predicted-taxid, "best" result for taxcheck
    for child in root[0]:
        if child.tag == 'predicted-taxid' and child.get('confidence') == 'HIGH':
            taxonomy['predicted'] = child.get('taxid')
            
    if taxonomy['predicted'] == taxonomy['best']:
        taxonomy.pop('best')
     
    return taxonomy, child_100

def run_DIAMOND(fasta: str, db: str, sensitivity: Literal["sensitive", "more-sensitive", "very-sensitive", "ultra-sensitive"] = "more-sensitive", coverage: float=95.0, threads: int=2, outdir: str=None, taxonlist: list[int]=None, outname: str='diamond'):
    outfile = Path(outdir, outname+'_'+sensitivity+'.tsv')
    logfile = Path(outdir, outname+'_'+sensitivity+'_logfile.txt')
    
    if taxonlist:
        taxonlist = [str(taxon) for taxon in taxonlist]
        taxonlist = ','.join(taxonlist)
        completed_blast = subprocess.run(
            ["diamond", "blastp", "-d", db, "-q", fasta, "--" + sensitivity, "--query-cover", str(coverage), "-p", str(threads), "-o", outfile, "--taxonlist", taxonlist,
                "--outfmt", str(6), "qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore"],
            shell=False, stderr=subprocess.PIPE, text=True)
    else:
        completed_blast = subprocess.run(
            ["diamond", "blastp", "-d", db, "-q", fasta, "--" + sensitivity, "--query-cover", str(coverage), "-p", str(threads), "-o", outfile,
                "--outfmt", str(6), "qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore"],
            shell=False, stderr=subprocess.PIPE, text=True)
    
    with open(logfile, "a") as f:
        f.write(completed_blast.stderr)

def read_DIAMOND_results(results: str, dir: str, cutoff: float=90.0) -> tuple[pd.core.frame.DataFrame, list[str], list[str]]:
    colnames = ['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore']
    
    data = pd.read_table(results, header=None, names=colnames, index_col='qseqid')
    data = data.sort_values(by=["qseqid", "pident", "length"], ascending=[True, False, False])

    # Deduplicate
    duplicates = data.loc[data.index.duplicated(), :]
    data = data.loc[~data.index.duplicated(), :]
    # Cutoff
    below_cutoff = data.loc[data['pident'] < cutoff]
    data = data.loc[data['pident'] >= cutoff]

    return data, duplicates, below_cutoff

def read_PGAP_results(results: str) -> list[str]:
    records = []
    for seq_record in SeqIO.parse(results, "fasta"):
        records.append(seq_record)

    return records

def get_mapping_table(mapping: list[str], diamond: str, pgap: str, dir: str, cutoff: float=90.0, database: Literal['nr', 'swissprot'] = 'nr') -> tuple[dict[str, int], list[str], list[str]]:
    pgap = read_PGAP_results(pgap)
    diamond, duplicates, below_cutoff = read_DIAMOND_results(diamond, dir, cutoff)
    
    no_diamond = []
    statistics = {
        'added': 0,
        'not_best_hits': 0,
        'not_found': 0
    }

    statistics['not_best_hits'] = len(duplicates)+len(below_cutoff)

    for record in tqdm(pgap):
        if record.name in diamond.index:
            new_id = diamond.loc[record.name, 'sseqid']
            if database == 'swissprot':
                new_id = new_id.split('|')[1]
            locus_tag = record.name
            locus_tag = 'G_' + locus_tag.replace('|','_').replace('.','_')
            mapping.append([locus_tag, new_id])
            statistics['added'] += 1
        else:
            no_diamond.append(record)
            statistics['not_found'] += 1
            
    SeqIO.write(no_diamond, str(Path(dir,f'DIAMOND_input_{database}.faa')), "fasta")

    return statistics, duplicates, below_cutoff

# TODO: integrate into the SPECIMEN reports class
class DiamondReport():
    def __init__(
        self,
        statistics: dict,
        duplicates: dict[str, list],
        below_cutoff: dict[str, list]
    ) -> None:
        self.statistics = statistics
        self.duplicates = duplicates
        self.below_cutoff = below_cutoff

    def visualise(self, color_palette: str = "YlGn") -> matplotlib.figure.Figure:
        try:
            cmap = matplotlib.colormaps[color_palette]
        except ValueError:
            warnings.warn('Unknown color palette, setting it to "YlGn"')
            cmap = matplotlib.colormaps["YlGn"]

        taxids = list(self.statistics.keys())

        fig, ax = plt.subplots()
        positions = np.arange(len(self.statistics[taxids[0]].keys()))

        if len(self.statistics.keys()) == 3:
            barWidth = 0.33
            bars1 = ax.bar(positions, self.statistics[taxids[0]].values(), label=list(self.statistics.keys())[0], width=0.3, color=cmap(0.25))
            ax.bar_label(bars1, self.statistics[taxids[0]].values())
            bars2 = ax.bar(positions+barWidth, self.statistics[taxids[1]].values(), label=list(self.statistics.keys())[1], width=0.3, color=cmap(0.5)) 
            ax.bar_label(bars2, self.statistics[taxids[1]].values())
            bars3 = ax.bar(positions+2*barWidth, self.statistics[taxids[2]].values(), label=list(self.statistics.keys())[2], width=0.3, color=cmap(0.75))
            ax.bar_label(bars3, self.statistics[taxids[2]].values())
            ax.set_ylabel('count', labelpad=10)
            plt.legend()
            plt.xticks(positions+barWidth, ['mapped', 'not best hits', 'not found'])
            plt.title('Report for the DIAMOND runs')

        elif len(self.statistics.keys()) == 2:
            barWidth = 0.45
            bars1 = ax.bar(positions, self.statistics[taxids[0]].values(), label=list(self.statistics.keys())[0], width=0.4, color=cmap(0.5)) 
            ax.bar_label(bars1, self.statistics[taxids[0]].values())
            bars2 = ax.bar(positions+barWidth, self.statistics[taxids[1]].values(), label=list(self.statistics.keys())[1], width=0.4, color=cmap(0.75))
            ax.bar_label(bars2, self.statistics[taxids[1]].values())
            ax.set_ylabel('count', labelpad=10)
            plt.legend()
            plt.xticks(positions+barWidth/2, ['mapped', 'not best hits', 'not found'])
            plt.title('Report for the DIAMOND runs')
        
        elif len(self.statistics.keys()) == 1:
            bars = ax.bar(self.statistics[taxids[0]].keys(), self.statistics[taxids[0]].values(), label=list(self.statistics.keys())[0], width=0.7, color=cmap(0.75))
            ax.bar_label(bars, self.statistics[taxids[0]].values())
            ax.set_ylabel('count', labelpad=10)
            plt.legend()
            plt.xticks(positions, ['mapped', 'not best hits', 'not found'])
            plt.title('Report for the DIAMOND run')
        
        return fig

    def save(self, dir: str, color_palette: str = 'YlGn') -> None:
        dir = str(Path(dir, 'DIAMOND_report'))
        Path(dir).mkdir(parents=True, exist_ok=False)

        fig = self.visualise(color_palette)
        fig.savefig(str(Path(dir, 'DIAMOND_visual.png')))
        
        for key in self.duplicates:
            pd.DataFrame(self.duplicates[key]).to_csv(str(Path(dir,f'duplicates_{key}.tsv')), sep="\t")
        for key in self.below_cutoff:
            pd.DataFrame(self.below_cutoff[key]).to_csv(str(Path(dir,f'below_cutoff_{key}.tsv')), sep="\t")

def adapt_config(cfg: dict[str, str], dir: str, modelname: str, refseq: str=None, genome: str=None) -> str:
    # needs to be added: gff, protein fasta, refseq accession, flag from_pgab
    config = validate_config(cfg['configpath'], cfg['next_step'])
    
    match cfg['next_step']:
        case 'cmpb':
            config['general']['gff'] = str(Path(dir,'pgap','output','annot.gff'))
            config['general']['protein_fasta'] = str(Path(dir,'pgap','output','annot_translated_cds.faa'))
            config['carveme']['refseq'] = refseq
            config['cm-polish']['is_lab_strain'] = True
        case 'hqtb':
            # TODO: integrate HQTB-pipeline -> fasta and gff are __USER__-parameters...
            # QUESTION: faa or gbff? both possible, test with faa
            # refineGEMS-Methode utility.io.mimic_genbank..., in case only gbff works
            if refseq:
                config['subject']['annotated_genome'] = cfg['fasta']
                config['subject']['gff'] = cfg['gff']
            else:
                config['subject']['annotated_genome'] = str(Path(dir,'pgap','output','annot_translated_cds.faa'))
                config['subject']['gff'] = str(Path(dir,'pgap','output','annot.gff'))
            config['subject']['full_sequence'] = genome
        case _:
            raise ValueError(f"Unknown input for pipeline: {cfg['next_step']}")
    config['general']['dir'] = dir
    config['general']['from_pgab'] = True
    
    new_configpath = str(Path(dir,f'{cfg['next_step']}_config.yaml'))
    with open(new_configpath, 'w') as f:
        yaml.dump(config, f)
    return new_configpath

# next methods