#!/usr/bin/env python3
"""Functions to run the workflow based on the Prokaryotic Genome Annotation Pipeline (PGAP). """

import datetime
import gzip
import logging
import matplotlib
import numpy as np
import os
import pandas as pd
import requests
import shutil
import subprocess
import tempfile
import time
import warnings
import yaml

from Bio import SeqIO
from importlib.resources import files
import matplotlib.pyplot as plt
from pathlib import Path
from tqdm.auto import tqdm
from typing import Union, Literal
import xml.etree.ElementTree as ET

from ..classes.reports import DIAMONDReport
from ..cmpb.workflow import run as run_cmpb
from ..hqtb.workflow import run as run_hqtb
from ..util.set_up import build_data_directories, validate_config

logger = logging.getLogger(__name__)

# TODO: improve variable names
# QUESTION: add logging from external methods to logging file?
# TODO: access points for command line
# TODO: PGAB wrapper
def run(configpath: Union[str, None] = None):
    """Run the PGAP-based (PGAB) workflow.

    Args:
        - configpath (Union[str,None]):
            The Path to a configuration file. If none given, prompts the user to 
            enter the needed information step by step.
            Defaults to None.
    """
    
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
            + str(datetime.date.today().year).removeprefix("20")
        )
    elif config["general"]["modelname"] is not None:
        modelname = config["general"]["modelname"]
    else:
        logger.info(
            "No values given for the standard name for a model. Default name will be used."
        )
        modelname = "model_" + str(datetime.date.today().year).removeprefix("20")

    # create log
    # ----------
    today = datetime.date.today().strftime("%Y%m%d")
    log_file = Path(dir, "logs", f"specimen_pgab_{str(today)}.log")
    handler = logging.handlers.RotatingFileHandler(
        log_file, mode="w", backupCount=10, encoding="utf-8", delay=0
    )
    handler.setFormatter(
        logging.Formatter(
            "{levelname} \t {name} \t {message}",
            style="{",
        )
    )
    logger.addHandler(handler)
    total_time_s = time.time()

    # run PGAP
    # --------
    logger.info('Running PGAP ...')
    step_start = time.time()
    
    # run_pgap(config['pgap'], dir)
    
    step_end = time.time()
    logger.info(f"\truntime: {step_end-step_start}s\n")

    match_100 = None
    match config['pgap']['tax-check']:
        case 'only':
            logger.info('The taxonomy check with PGAP is finished. The results can be viewed in XML- or txt-format.')
            result = parse_tax_check(str(Path(dir,'pgap','output','ani-tax-report.xml')), type='only')
            logger.info('As a further step, running the HQTB pipeline with one of these strains as a template is recommended:')
            logger.info(result.to_string(index=False))
            return
        case 'continue':
            taxonomy, match_100 = parse_tax_check(str(Path(dir,'pgap','output','ani-tax-report.xml')), config['pgap']['tax_cutoff'])
            # transform dict into a list for DIAMOND
            if len(taxonomy) > 1:
                taxid = [taxonomy['best'], taxonomy['predicted']]
            else: taxid = [taxonomy['predicted']]
            
            # TODO: default cutoff?
            if taxid:
                if config['pgap']['diamond']=='run':
                    logger.info(f'The taxid(s) which will be used for the taxonomy specific DIAMOND run are {taxid}.')
            else:
                logger.warning('There was no taxid found above the cutoff. You will have to decide on a taxid manually.')
                return
        case 'none':
            if config['pgap']['taxid']:
                taxid = str(config['pgap']['taxid']).split(',')
                taxid = [t.strip() for t in taxid]
            else: taxid = None
        case _:
            raise ValueError(f"Unknown input for taxonomy check: {config['pgap']['tax-check']}")

    # option to skip DIAMOND if a 100% match is found
    if match_100 and config['pgap']['diamond']!='run':
        # TODO: maybe revise the structure?
        logger.info(f'DIAMOND will be skipped because of a found organism with 100% identity. As a next step, the {config['polish']['next_step']} workflow will be started.')

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
            match config['polish']['next_step'].lower():
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
                        logger.warning(f"There was no protein FASTA file found for {refseq} the automated way. To run HQTB, you will have to specify the files yourself.")
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
                        logger.warning(f"There was no GFF file found for {refseq} the automated way. To run HQTB, you will have to specify the files yourself.")
                        return
                    polish_config = adapt_config(config['polish'], dir, modelname, refseq, genome=config['pgap']['generic']['fasta']['location'])
                    run_hqtb(polish_config)
                case 'stop':
                    logger.warning('There was no following polishing step defined. This workflow ends here.')
                case _:
                    raise ValueError(f"Unknown input for next polishing step: {config['polish']['next_step']}")
            return
        else: 
            logger.warning('There was no matching RefSeq-ID found and CarveMe could not be performed.')
            if config['pgap']['diamond']=='maybe':
                logger.info(f'Because of this, DIAMOND will run with the following taxonomy id(s): {taxid}')
            else: return

    # data structures for nr BLAST
    fasta_nr = str(Path(dir,'pgap','output','annot_translated_cds.faa'))
    mapping_nr = []
    # information for the DIAMONDReport
    duplicates_nr = {}
    below_cutoff_nr = {}
    statistics_nr = {}
    
    # data structures for SwissProt BLAST
    fasta_sp = str(Path(dir,'pgap','output','annot_translated_cds.faa'))
    mapping_sp = []
    # information for the DIAMONDReport
    duplicates_sp = {}
    below_cutoff_sp = {}
    statistics_sp = {}
    
    # data structures for USER BLAST
    fasta_ub = str(Path(dir,'pgap','output','annot_translated_cds.faa'))
    mapping_ub = []
    # information for the DIAMONDReport
    duplicates_ub = {}
    below_cutoff_ub = {}
    statistics_ub = {}

    if taxid is not None:
        for id in taxid:
            run_id = []
            run_id.append(id)
            
            # nr database
            if config['nr_blast']['run']:
                logger.info(f'Running DIAMOND run with taxid {id} and nr database ...')
                step_start = time.time()
                
                run_DIAMOND(fasta_nr, config['nr_blast']['db'], 
                            config['nr_blast']['sensitivity'], config['nr_blast']['coverage'], config['nr_blast']['threads'], 
                            str(Path(dir,'nr_blast')), run_id, f'nr_{id}')
                stats_tax, duplicates_tax, cutoff_tax = get_mapping_table(mapping_nr, str(Path(dir,'nr_blast',f'nr_{id}_{config['nr_blast']['sensitivity']}.tsv')), fasta_nr, dir, config['nr_blast']['pid'], 'nr')
                statistics_nr[id] = stats_tax
                duplicates_nr[id] = duplicates_tax
                below_cutoff_nr[id] = cutoff_tax
                
                step_end = time.time()
                logger.info(f"\truntime: {step_end-step_start}s\n")

            fasta_nr = str(Path(dir,'remaining_proteins_nr.faa'))
            
            # SwissProt
            if config['swissprot_blast']['run']:
                logger.info(f'Running DIAMOND run with taxid {id} and SwissProt database ...')
                step_start = time.time()
                
                run_DIAMOND(fasta_sp, config['swissprot_blast']['db'], 
                            config['swissprot_blast']['sensitivity'], config['swissprot_blast']['coverage'], config['swissprot_blast']['threads'], 
                            str(Path(dir,'swissprot_blast')), run_id, f'swissprot_{id}')
                stats_tax, duplicates_tax, cutoff_tax = get_mapping_table(mapping_sp, str(Path(dir,'swissprot_blast',f'swissprot_{id}_{config['swissprot_blast']['sensitivity']}.tsv')), fasta_sp, dir, config['swissprot_blast']['pid'], 'swissprot')
                statistics_sp[id] = stats_tax
                duplicates_sp[id] = duplicates_tax
                below_cutoff_sp[id] = cutoff_tax
                
                step_end = time.time()
                logger.info(f"\truntime: {step_end-step_start}s\n")
            
            fasta_sp = str(Path(dir,'remaining_proteins_swissprot.faa'))
            
            # USER database
            if config['user_blast']['run']:
                logger.info(f'Running DIAMOND run with taxid {id} and USER database ...')
                step_start = time.time()
                
                run_DIAMOND(fasta_ub, config['user_blast']['db'], 
                            config['user_blast']['sensitivity'], config['user_blast']['coverage'], config['user_blast']['threads'], 
                            str(Path(dir,'user_blast')), run_id, f'user_{id}')
                stats_tax, duplicates_tax, cutoff_tax = get_mapping_table(mapping_ub, str(Path(dir,'user_blast',f'user_{id}_{config['user_blast']['sensitivity']}.tsv')), fasta_ub, dir, config['user_blast']['pid'], 'user')
                statistics_ub[id] = stats_tax
                duplicates_ub[id] = duplicates_tax
                below_cutoff_ub[id] = cutoff_tax
                
                step_end = time.time()
                logger.info(f"\truntime: {step_end-step_start}s\n")
            
            fasta_ub = str(Path(dir,'remaining_proteins_user.faa'))

    # run DIAMOND without taxonomy sensitivity
    if config['nr_blast']['run']:
        logger.info(' Running DIAMOND run without taxonomy restriction and nr database ...')
        step_start = time.time()
        
        run_DIAMOND(fasta_nr, config['nr_blast']['db'], 
                    config['nr_blast']['sensitivity'], config['nr_blast']['coverage'], config['nr_blast']['threads'], 
                    str(Path(dir,'nr_blast')), outname='nr_not_taxrestricted')
        
        step_end = time.time()
        logger.info(f"\truntime: {step_end-step_start}s\n")

        # Expand mapping table
        stats_notax, duplicates_notax, cutoff_notax = get_mapping_table(mapping_nr, str(Path(dir,'nr_blast','nr_not_taxrestricted_'+config['nr_blast']['sensitivity']+'.tsv')), fasta_nr, dir, config['nr_blast']['pid'], 'nr')
        statistics_nr['no_tax_restriction'] = stats_notax
        duplicates_nr['no_tax_restriction'] = duplicates_notax
        below_cutoff_nr['no_tax_restriction'] = cutoff_notax
        
        # DIAMOND report
        report = DIAMONDReport(statistics_nr, duplicates_nr, below_cutoff_nr)
        report.save(str(Path(dir,'nr_blast')))
        
    if config['swissprot_blast']['run']:
        logger.info('Running DIAMOND run without taxonomy restriction and SwissProt database ...')
        step_start = time.time()
        
        run_DIAMOND(fasta_sp, config['swissprot_blast']['db'], 
                    config['swissprot_blast']['sensitivity'], config['swissprot_blast']['coverage'], config['swissprot_blast']['threads'], 
                    str(Path(dir,'swissprot_blast')), outname='swissprot_not_taxrestricted')
        
        step_end = time.time()
        logger.info(f"\truntime: {step_end-step_start}s\n")

        # Expand mapping table
        stats_notax, duplicates_notax, cutoff_notax = get_mapping_table(mapping_sp, str(Path(dir,'swissprot_blast','swissprot_not_taxrestricted_'+config['swissprot_blast']['sensitivity']+'.tsv')), fasta_sp, dir, config['swissprot_blast']['pid'], 'swissprot')
        statistics_sp['no_tax_restriction'] = stats_notax
        duplicates_sp['no_tax_restriction'] = duplicates_notax
        below_cutoff_sp['no_tax_restriction'] = cutoff_notax
        
        # DIAMOND report
        report = DIAMONDReport(statistics_sp, duplicates_sp, below_cutoff_sp)
        report.save(str(Path(dir,'swissprot_blast')))
        
    if config['user_blast']['run']:
        logger.info('Running DIAMOND run without taxonomy restriction and USER database ...')
        step_start = time.time()
        
        run_DIAMOND(fasta_ub, config['user_blast']['db'], 
                    config['user_blast']['sensitivity'], config['user_blast']['coverage'], config['user_blast']['threads'], 
                    str(Path(dir,'user_blast')), outname='user_not_taxrestricted')
        
        step_end = time.time()
        logger.info(f"\truntime: {step_end-step_start}s\n")

        # Expand mapping table
        stats_notax, duplicates_notax, cutoff_notax = get_mapping_table(mapping_ub, str(Path(dir,'user_blast','user_not_taxrestricted_'+config['user_blast']['sensitivity']+'.tsv')), fasta_ub, dir, config['user_blast']['pid'], 'user')
        statistics_ub['no_tax_restriction'] = stats_notax
        duplicates_ub['no_tax_restriction'] = duplicates_notax
        below_cutoff_ub['no_tax_restriction'] = cutoff_notax
        
        # DIAMOND report
        report = DIAMONDReport(statistics_ub, duplicates_ub, below_cutoff_ub)
        report.save(str(Path(dir,'user_blast')))
        
    # save mapping tables (even if empty)
    mapping_nr = pd.DataFrame(mapping_nr, columns=['model_id','NCBI'])
    mapping_nr.to_csv(str(Path(dir, 'mapping_table_nr.csv')), header=True, index=False, sep='\t')
    # QUESTION: UniProt instead of DECLASSIFIED?
    mapping_sp = pd.DataFrame(mapping_sp, columns=['model_id','UNIPROT'])
    mapping_sp.to_csv(str(Path(dir, 'mapping_table_swissprot.csv')), header=True, index=False, sep='\t')
    # USER database: has to be handled manually
    logging.info('The results from user specific BLAST have to be manually added to the model.')
    mapping_ub = pd.DataFrame(mapping_sp, columns=['model_id','UNCLASSIFIED'])
    mapping_ub.to_csv(str(Path(dir, 'mapping_table_user.csv')), header=True, index=False, sep='\t')

    logging.info('The remaining proteins files are proteins that could not be blasted successfully.')
    
    polish_config = adapt_config(config['polish'], dir, modelname, genome=config['pgap']['generic']['fasta']['location'])
    match config['polish']['next_step'].lower():
        case 'cmpb':
            logger.info(f'Running next polishing step {config['polish']['next_step']} ...')
            step_start = time.time()
            
            run_cmpb(polish_config)
            # print('work in progress...')
        case 'hqtb':
            logger.info(f'Running next polishing step {config['polish']['next_step']} ...')
            step_start = time.time()
            
            # TODO: automated way to find a template model?
            run_hqtb(polish_config)
            # print('work in progress...')
        case 'stop':
            logger.warning('There was no following polishing step defined. This workflow ends here.')
        case _:
            raise ValueError(f"Unknown input for next polishing step: {config['polish']['next_step']}")
        
    step_end = time.time()
    logger.info(f"\truntime: {step_end-step_start}s\n")
    
    total_time_e = time.time()
    logger.info(f'Total run time: {total_time_e - total_time_s}')

    # QUESTION: What exactly belongs in the header of the mapping table? 'NCBI' for protein_tag and 'UNIPROT' for UniProt?
    
# TODO
def save_pgab_user_input() -> dict[str, str]:
    """method adapted from SPECIMEN.util.set_up
    !under construction!
    """

def run_pgap(config: dict, dir: str):
    """Run PGAP on the command line. Saves the results in a folder called 'output'.

    Args:
        - config (dict):
            Dictionary with input parameters for PGAP.
            Needs to have the keys 'generic' and 'metadata'.
        - dir (str):
            Path to the directory where the results are saved to.

    Raises:
        ValueError: Unknown input for taxonomy check parameter.
    """
    
    logfile = str(Path(dir,'logs','log_pgap.txt'))
    dir_pgap = str(Path(dir,'pgap'))

    prepare_pgap_input(config['generic'], config['metadata'], dir_pgap)

    match config['tax-check']:
        # BUG: does not work without --ignore-all-errors
        case 'only': # user searches an organism, no PGAP
            completed_pgap = subprocess.run(['./pgap.py', '-r', '--taxcheck-only', '--no-self-update', '--ignore-all-errors', '-o', str(Path(dir_pgap, 'output')), str(Path(dir_pgap, 'input.yaml'))],
                shell=False, stderr=subprocess.PIPE, text=True)
        case 'continue': # user does taxcheck and PGAP
            completed_pgap = subprocess.run(['./pgap.py', '-r', '--taxcheck', '--no-self-update', '--ignore-all-errors', '-o', str(Path(dir_pgap, 'output')), str(Path(dir_pgap, 'input.yaml'))],
                shell=False, stderr=subprocess.PIPE, text=True)
        case 'none': # user already has an organism
            completed_pgap = subprocess.run(['./pgap.py', '-r', '--no-self-update', '--ignore-all-errors', '-o', str(Path(dir_pgap, 'output')), str(Path(dir_pgap, 'input.yaml'))],
                shell=False, stderr=subprocess.PIPE, text=True)
        case _:
            raise ValueError(f"Unknown input for taxonomy check: {config['pgap']['tax-check']}")
    with open(logfile, "a") as f:
        f.write(completed_pgap.stderr)

def prepare_pgap_input(generic: dict, metadata: dict, dir: str):
    """Function to convert PGAB config to input files for PGAP.

    Args:
        - generic (dict): 
            Dictionary with the path to the genomic FASTA and
            the submol.yaml file.
        - metadata (dict): 
            Dictionary with metadata for PGAP.
        - dir (str): 
            Path to the directory where the files are saved to.
    """
    
    subprocess.run(['cp', generic['fasta']['location'], dir])
    
    if 'position' in metadata.keys:
        metadata['location'] = metadata['position']
        metadata.pop('position')
    else: logger.warning('KeyError for \'location\' in metadata dictionary')

    with open(Path(dir,'submol.yaml'), 'w') as f:
        yaml.dump(metadata, f)

    generic.update({'submol': {'class': 'File', 'location': str(Path(dir,'submol.yaml'))}})
    generic['fasta']['location'] = generic['fasta']['location'].split('/')[-1]

    with open(Path(dir,'input.yaml'), 'w') as f:
        yaml.dump(generic, f)

def parse_tax_check(file: str, cutoff: float=90.0, type: Literal['only','continue'] = 'continue') -> tuple[dict[str, str], dict[str, str]|None]:
    """Parses the tax report XML file from PGAP.

    Args:
        - file (str):
            Path to the tax report from PGAP.
        - cutoff (float, optional): 
            Cutoff for organism similarity.
            Defaults to 90.0.
        - type (Literal["only","continue"], optional): 
            Describes if only the taxonomy check ran.
            Defaults to 'continue'.

    Returns:
        dict[str, str]:
            Dictionary with the taxonomy results above the threshold.
        dict[str, str]|None:
            Dictionary of a 100% match if one was found.
    """
    
    taxonomy = {}
    results = []
    child_100 = {}

    data = ET.parse(file)
    root = data.getroot()

    for child in root[1]:
        if child.tag == 'subject' :
            # Search for a 100% match
            if float(child.get('ANI')) == 100:
                logger.info(f"An organism with a 100% identity was found. This organism is {child.get('org-name')}.")
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

def run_DIAMOND(fasta: str, db: str, sensitivity: str = "more-sensitive", # QUESTION: alle Optionen in Literal[]?
                coverage: float=95.0, threads: int=2, outdir: str=None, taxonlist: list[int]=None, outname: str='diamond'):
    """Run DIAMOND BLASTp on the command line.

    Args:
        - fasta (str): 
            Path to the FASTA input file.
        - db (str):
            Path to the dmnd database.
        - sensitivity (optional): 
            Sensitivity option for DIAMOND.
            Defaults to "more-sensitive".
        - coverage (float, optional): 
            Coverage for DIAMOND. Defaults to 95.0.
        - threads (int, optional): 
            Defaults to 2.
        - outdir (str, optional):
            Path to directory to save the DIAMOND result file..
        - taxonlist (list[int], optional): 
            List of taxonomy IDs for taxonomy specific BLAST runs.
            Defaults to None.
        - outname (str, optional): 
            File name for the output and logging files. Defaults to 'diamond'.
    """
    
    outfile = Path(outdir, outname+'_'+sensitivity+'.tsv')
    dir_log = str(Path(Path(outdir).parents[0],'logs'))
    logfile = Path(dir_log, 'log_'+outname+'_'+sensitivity+'.txt')
    
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

def read_DIAMOND_results(results: str, cutoff: float=90.0) -> tuple[pd.core.frame.DataFrame, list[str], list[str]]:
    """Function to read in the DIAMOND results and sort out
    duplicates and hits below the cutoff.

    Args:
        - results (str):
            Path to the DIAMOND tsv file.
        - cutoff (float, optional): 
            Cutoff for the percentage identity. Defaults to 90.0.

    Returns:
        pd.core.frame.DataFrame:
            Dataframe with the found BLAST hits.
        list[str]:
            List with the found duplicates.
        list[str]:
            List with the hits below the cutoff.
    """
    
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
    """Function to read in the PGAP FASTA file.

    Args:
        - results (str):
            Path to the FASTA file.

    Returns:
        list[str]:
            List with the annotated proteins from PGAP.
    """
    
    records = []
    for seq_record in SeqIO.parse(results, "fasta"):
        records.append(seq_record)

    return records

def get_mapping_table(mapping: list[str], diamond: str, pgap: str, dir: str, cutoff: float=90.0, database: Literal['nr', 'swissprot'] = 'nr') -> tuple[dict[str, int], list[str], list[str]]:
    """Function to combine DIAMOND and PGAP results as an ID mapping table.

    Args:
        - mapping (list[str]):
            (Empty) List as the mapping table.
        - diamond (str):
            Path to the DIAMOND tsv file.
        - pgap (str):
            Path to the PGAP FASTA file.
        - dir (str):
            Path to the directory where the remaining proteins are saved to
            which where not found with DIAMOND.
        - cutoff (float, optional): 
            Cutoff for the percentage identity. Defaults to 90.0.
        - database (Literal["nr", "swissprot"], optional): 
            Database with which the BLASTp run was performed.
            Defaults to 'nr'.

    Returns:
        dict[str, int]:
            Dictionary with 'mapped', 'not best hits' and 'not found' as the keys
            and the respective counts as values.
        list[str]:
            List of the duplicates found with DIAMOND.
        list[str]:
            List of DIAMOND hits below the cutoff.
    """
    
    pgap = read_PGAP_results(pgap)
    diamond, duplicates, below_cutoff = read_DIAMOND_results(diamond, cutoff)
    
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
            
    SeqIO.write(no_diamond, str(Path(dir,f'remaining_proteins_{database}.faa')), "fasta")

    return statistics, duplicates, below_cutoff

def adapt_config(cfg: dict[str, str], dir: str, modelname: str, refseq: str=None, genome: str=None) -> str:
    """Function to adapt a CMPB or HQTB config with new parameters.

    Args:
        cfg (dict[str, str]): 
            Dictionary with information for the file to be adapted.
        dir (str):
            Path to the directory to save the config to.
        modelname (str):
            String to be used as the model name in the next steps.
        refseq (str, optional): 
            RefSeq accession if known.
            Defaults to None.
        genome (str, optional): 
            Path to the genomic FASTA file. 
            Defaults to None.

    Raises:
        ValueError: Unknown input for pipeline.

    Returns:
        str:
            Path to the new config.
    """
    
    if cfg['next_step'].lower()=='stop':
        return
    
    config = validate_config(cfg['configpath'], cfg['next_step'])
    
    match cfg['next_step'].lower():
        case 'cmpb':
            config['general']['gff'] = str(Path(dir,'pgap','output','annot.gff'))
            config['general']['protein_fasta'] = str(Path(dir,'pgap','output','annot_translated_cds.faa'))
            config['carveme']['refseq'] = refseq
            config['cm-polish']['is_lab_strain'] = True
        case 'hqtb':
            # QUESTION: pseudo paths for the user input? __USER__
            # yaml.load instead of pseudo paths?
            if refseq:
                config['subject']['annotated_genome'] = cfg['fasta']
                config['subject']['gff'] = cfg['gff']
            else:
                config['subject']['annotated_genome'] = str(Path(dir,'pgap','output','annot_translated_cds.faa'))
                config['subject']['gff'] = str(Path(dir,'pgap','output','annot.gff'))
            config['subject']['full_sequence'] = genome
        case 'stop':
            logger.warning('There was no following polishing step and no config defined . This workflow ends here.')
            return
        case _:
            raise ValueError(f"Unknown input for pipeline: {cfg['next_step']}")
    config['general']['dir'] = dir
    config['general']['from_pgab'] = True
    
    new_configpath = str(Path(dir,f'{cfg['next_step']}_config.yaml'))
    with open(new_configpath, 'w') as f:
        yaml.dump(config, f)
    return new_configpath

# next methods