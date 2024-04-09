"""Part one of the third step of the pipeline: refinement - extension.

Extents the model by mapping missing genes to reactions via NCBI, KEGG,
MetaNetX and BiGG.
"""
__author__ = 'Carolin Brune'
################################################################################
# requirements
################################################################################

from pathlib import Path
from Bio import SeqIO

import pandas as pd
import re
import subprocess
import sys
import time
import urllib.error

from tqdm import tqdm
from tqdm.auto import tqdm
tqdm.pandas()
import cobra

from typing import Literal

from Bio.KEGG import REST
from Bio.KEGG import Enzyme, Compound

from refinegems.utility.io import kegg_reaction_parser, load_a_table_from_database
from refinegems.utility.entities import create_random_id, get_reaction_annotation_dict, match_id_to_namespace
from refinegems.analysis.investigate import run_memote

# further required programs:
#        - DIAMOND, tested with version 0.9.14 (works only for certain sensitivity mode)
#                   tested with version 2.1.8 (works for all sensitivity modes for that version)

################################################################################
# functions
################################################################################

# ----------------------
# identify missing genes
# ----------------------
def identify_missing_genes(gene_path:str, model:cobra.Model, 
                           id:str, fasta_path:str, dir:str) -> pd.DataFrame:
    """Identify missing genes in a model based on the draft model and a list of genes.
    Save the missing genes' sequences in a new FASTA file.


    Args:
        - gene_path (str): Path to the csv file containing the list of genes.
        - model (cobra.Model): The draft model.
        - id (str): Name of the column in gene_path, that equal the ids in the model.
        - fasta_path (str): Path to the FASTA file containing all genes for the genome.
        - dir (str): Path to the directory to save the new FASTA in.

    Returns:
        pd.DataFrame:  A table with the missing genes' information.
    """

    # get list of all genes
    gene_list_all = pd.read_csv(gene_path)
    # filter list by genes in model
    gene_list_model = [g.id for g in model.genes]
    gene_list_missing = gene_list_all[gene_list_all[id].isin(gene_list_model) == False]
    # extract wanted sequences from the FASTA and save new FASTA
    fasta = SeqIO.parse(fasta_path, 'fasta')
    SeqIO.write((seq for seq in fasta if seq.id in gene_list_missing[id].tolist()), Path(dir,'missing_genes.fasta'), 'fasta')

    return gene_list_missing


def find_best_diamond_hits(file:str, pid:float) -> pd.DataFrame:
    """Identify the best hits from the DIAMOND run.

    Args:
        - file (str): File name (path) of the "missing_genes" file.
        - pid (float): Threshold value for the PID. 
            Should be between 0.0 and 100.0 (given in percentage)

    Returns:
        pd.DataFrame: Table of the locus tags of the query and the NCBI accession version numbers of the RefSeqs.
    """

    # load diamond results
    diamond_results = pd.read_csv(file, sep='\t')
    diamond_results.columns = ['query_ID', 'subject_ID', 'PID', 'align_len', 'no_mismatch', 'no_gapopen', 'query_start', 'query_end', 'subject_start', 'subject_end','E-value','bitscore']
    # filter by PID
    diamond_results = diamond_results[diamond_results['PID']>=pid]
    # get best hit for each gene in list
    # ...............................................
    # @TODO
    # NOTE: best hit = FIRST best hit for a locus_tag
    #       mainly a runtime optimisation
    #       however can lead to information loss
    # ...............................................
    best_hits = diamond_results.loc[diamond_results.groupby('query_ID')['PID'].idxmax()].iloc[:,[0,1]]
    best_hits.columns = ['locus_tag','ncbi_accession_version']

    return best_hits

# -----------------------------
# map to additional information
# -----------------------------

# mapping zo NCBI
# ---------------
def map_to_NCBI_efetch_row(row:pd.Series) -> pd.Series:
    """Map a single entry from the table in get_ncbi_info() to NCBI using EntrezDirect.

    Args:
        - row (pd.Series): A single row of the table.

    Returns:
        pd.Series: The edited row.
    """

    #  EC number
    # fetch_ec = subprocess.run(['efetch','-db','protein','-id',row["ncbi_accession_version"],'-format','gpc','|','xtract','-insd','Protein','EC_number'], stdout=subprocess.PIPE, shell=False, encoding='utf-8')
    fetch_ec = subprocess.run(['efetch','-db','protein','-id',row["ncbi_accession_version"],'-format','gpc'], check=True, capture_output=True)
    fetch_ec = subprocess.run(['xtract','-insd','Protein','EC_number'], input=fetch_ec.stdout, capture_output=True)
    fetch_ec = fetch_ec.stdout.decode('utf-8')
    if fetch_ec:
        row['EC number'] = fetch_ec
    else:
        row['EC number'] = '-'

    # locus_tag_ref old_locus_tag GeneID
    fetch_rest = subprocess.run(['esearch','-db','gene','-query',row["ncbi_accession_version"]], check=True, capture_output=True)
    fetch_rest = subprocess.run(['efetch','-format','docsum'], input=fetch_rest.stdout, capture_output=True)
    fetch_rest = subprocess.run(['xtract','-pattern','DocumentSummary','-element','Id,OtherAliases'], input=fetch_rest.stdout, capture_output=True)
    fetch_rest = fetch_rest.stdout.decode('utf-8')

    if fetch_rest:
        # since there is no way to determine, which entry is the correct one,
        # if multiple entries are returned, always take the first one.
        # Assuming, if multiple are returned, they represent the same gene in different subjects.
        fetch_rest = fetch_rest.split('\n')[0]
        id, aliases = fetch_rest.split('\t')
        row['GeneID'] = id
        aliases = [a.strip for a in aliases.split(',')]
        row['locus_tag_ref'] = aliases[0]
        if len(aliases) > 1:
            row['old_locus_tag'] = aliases[1]
            # further names get ignored
            # .........................................
            # possible problem:
            # assumes, that aliases 1 and 2 are always
            # locus_tag_ref and old_locus_tag
            # .........................................

    return row


# @TODO
def get_ncbi_info(table:pd.DataFrame, ncbi_map_file:str=None) -> pd.DataFrame:
    """Retrieve information from NCBI via mapping to precomputed files.
       If no mapping exists, try via Entrez - implemented, but does not work well
       -> problems with connecting successfully to NCBI.
       @TODO if time, rewrite using the Bioppython package, see if it works better

    Args:
        - table (pd.DataFrame):  A table containing the output of "find_best_diamond_hits".
            Must have the column 'ncbi_accession_version'.
        - ncbi_map_file (str, optional): A csv file containing the ncbi mapping. 
            Defaults to None.

    Returns:
        pd.DataFrame:  A table with the added information.
    """

    # case 1: if precomputed mapping was given
    if ncbi_map_file:
        print('\t\tRetrieving mapping from precomputed file.')
        # load pre-compiles ncbi mapping
        ncbi_mapping = pd.read_csv(ncbi_map_file)

        # merge with input info table:
        # locus_tag ncbi_accession_version locus_tag_ref old_locus_tag GeneID EC number
        # ......................................................................
        # @TODO // @NOTE
        #    using first() takes the first mapped accession version organism
        #    instead of using all
        #    while testing, all matches returned same EC number, which is why
        #    this is used to get a better runtime
        #    however, if the locus tag play a more significant role,
        #    it might be better to remove the '.first().reset_index()'
        # ......................................................................
        table = pd.merge(table, ncbi_mapping, on='ncbi_accession_version').fillna('-').groupby('locus_tag').first().reset_index()

    # case 2: no mapping exists
    else:
        print('\t\tNo precomputed mapping. Retrieving information directly from NCBI.\n\t\tThis may take a while.')
        table['locus_tag_ref'] = pd.Series(dtype='str')
        table['old_locus_tag'] = pd.Series(dtype='str')
        table['GeneID'] = pd.Series(dtype='str')
        table['EC number'] = pd.Series(dtype='str')
        table = table.progress_apply(lambda row: map_to_NCBI_efetch_row(row), axis=1)

    return table

# mapping to KEGG
# ---------------
# @TODO
def map_to_KEGG_row(row:pd.Series, loci_list:list[str]=None) -> pd.Series:
    """Map a single entry of the table in map_to_KEGG() to KEGG.

    Args:
        - row (pd.Series): A single row of the table.
        - loci_list (list[str], optional): List of starting letters of locus tags belong to organisms, that are listed in KEGG
            and the corresponding organism three-letter code. 
            Defaults to None (= search all entries).

    Returns:
        pd.Series: The mapped row.
    """

    # --------------------------------
    # decide search based on loci list
    # --------------------------------
    # @TODO if no list is given, all locus tags will be searched
    # if no list is given, all locus tags will be searched
    search_locus = False
    # if a list of locus tags in KEGG is given
    if loci_list:
        # use it to search KEGG more efficiently
        if any(x in row['locus_tag_ref'] for x in loci_list):
            # search only if loci might indeed be in KEGG
            search_locus = True
            #@TODO this assumes that the locus tag is of the format XXX_YYY
            # can this be written more general?
            organism_abbrv = loci_list[row['locus_tag_ref'].split('_')[0]]+':'
        else:
            # if locus not in list, do not search
            search_locus = False

    # check, if ...
    # ----------------------------------------------
    # 1) ... organism (= locus tag start) is in KEGG
    # ----------------------------------------------
    if search_locus:
        try:
            # search locus tag in KEGG
            request = REST.kegg_get(organism_abbrv+row['locus_tag_ref'])
            for line in request.read().split('\n'):
                if line.startswith('ORTHOLOGY'):
                    if re.search(r'\[EC:(.*?)\]', line):
                        # @TODO
                        # will this result in problems with multiple EC numbers ?
                        ec = re.search(r'\[EC:(.*?)\]',line).group(1)
                        row['EC number'] = ec.strip()
                        row['KEGG.notes'] = 'updated EC with KEGG: True'
        except urllib.error.HTTPError:
            # since locus tag was not found, try ...
            # 1. : the old one (if it exists)
            found = False
            not_found_log = F'HTTPError: {organism_abbrv+row["locus_tag_ref"]}'
            if row['old_locus_tag'] != '-' and row['old_locus_tag'] != row['locus_tag_ref']:
                try:
                    request = REST.kegg_get(organism_abbrv+row['old_locus_tag'])
                    found = True
                    for line in request.read().split('\n'):
                        if line.startswith('ORTHOLOGY'):
                            if re.search(r'\[EC:(.*?)\]', line):
                                ec = re.search(r'\[EC:(.*?)\]',line).group(1)
                                row['EC number'] = ec.strip()
                                row['KEGG.notes'] = 'updated EC with KEGG: True'
                except urllib.error.HTTPError:
                    # if there is still a 'not found' error:
                    not_found_log += F', {organism_abbrv+row["old_locus_tag"]}'
                except ConnectionResetError:
                    print(F'ConnectionResetError: {organism_abbrv+row["old_locus_tag"]}')
                except urllib.error.URLError:
                    print(F'URLError: {organism_abbrv+row["old_locus_tag"]}')
            # 2. : check if RS after _, remove it and try again
            if '_RS' in row['locus_tag_ref'] and not found and row['old_locus_tag'] != row['locus_tag_ref'].replace('_RS','_'):
                try:
                    edited_tag = row['locus_tag_ref'].replace('_RS','_')
                    request = REST.kegg_get(organism_abbrv+edited_tag)
                    for line in request.read().split('\n'):
                        if line.startswith('ORTHOLOGY'):
                            if re.search(r'\[EC:(.*?)\]', line):
                                ec = re.search(r'\[EC:(.*?)\]',line).group(1).strip()
                                row['EC number'] = ec.strip()
                                row['KEGG.notes'] = 'updated EC with KEGG: True'
                    found = True
                except urllib.error.HTTPError:
                    # if there is still a 'not found' error:
                    not_found_log += F', {organism_abbrv+edited_tag}'
                except ConnectionResetError:
                    print(F'ConnectionResetError: {organism_abbrv+edited_tag}')
                except urllib.error.URLError:
                    print(F'URLError: {organism_abbrv+edited_tag}')
            if not found:
                print(not_found_log)
        except ConnectionResetError:
            print(F'ConnectionResetError: {organism_abbrv+row["locus_tag_ref"]}')
        except urllib.error.URLError:
            print(F'URLError: {organism_abbrv+row["locus_tag_ref"]}')

    # -------------------------------------------------
    # 2) ... EC number is in KEGG and leads to reaction
    # -------------------------------------------------
    if not '-' in row['EC number']:
        # .......................................
        # @TODO
        #    how to handle multiple EC numbers ?
        #    idea: split row but how ?
        # .......................................
        if ' ' in row['EC number']:
            row['EC number'] = row['EC number'].split(' ')[0]
            print(F'NOTE: Dropped EC numbers {row["EC number"].split(" ")[1:]} for {row["locus_tag"]}.')
        try:
            request = REST.kegg_get(F'ec:{row["EC number"]}')
            record = Enzyme.read(request)

            if len(record.reaction) == 0 or record.reaction == None:
                # no reaction associated
                return row

            elif len(record.reaction) == 1:
                # exactly one reaction (block) ...
                if re.search(r'\[RN:(.*?)\]', record.reaction[0]):
                    # ... with an ID in KEGG
                    row['KEGG.reaction'] = re.search(r'\[RN:(.*?)\]', record.reaction[0]).group(1)

            else:
                # multiple reaction(s) blocks found, no clear 'answer', which is the correct one
                # ......................
                # @TODO
                # current solution: skip
                # ......................
                print(F'Warning: multiple entries found for {row["EC number"]}, assignment skipped')
        except urllib.error.HTTPError:
            print(F'HTTPError: {row["EC number"]}')
        except ConnectionResetError:
            print(F'ConnectionResetError: {row["EC number"]}')
        except urllib.error.URLError:
            print(F'URLError: {row["EC number"]}')

    return row


def map_to_KEGG(gene_table:pd.DataFrame,working_dir:str,
                manual_dir:str,genome_summary:str=None) -> Path:
    """Map entries of the table from get_ncbi_info() to KEGG.
    Get KEGG.reaction ID from locus tag, if available, otherwise try EC number.
    Entries of the finalised table are saved depending whether they were assigned
    a KEGG.reaction ID (further use) or not (possible starting point for manual curation).

    Args:
        - gene_table (pd.DataFrame): The input table containing the output of get_ncbi_info().
        - working_dir (str): Path to the working directory - place to save genes with a KEGG.reaction ID.
        - manual_dir (str): Path to a directory to save the genes without KEGG assignment for possible manual curation.
        - genome_summary (str, optional):  Path to a information file for the reference genomes to check if they are in KEGG. 
            Defaults to None.

    Returns:
        Path: The path to the file containing the genes with KEGG assignment.
    """

    if genome_summary:
        # identify locus tags that belong to genomes,
        # that are listed in KEGG
        genome_info = pd.read_csv(genome_summary)
        # loci_in_KEGG = genome_info.loc[genome_info['KEGG.organism'].notna()]['locus_tag (start)'].tolist()
        loci_in_KEGG = genome_info.loc[genome_info['KEGG.organism'].notna()][['locus_tag (start)', 'KEGG.organism']].set_index('locus_tag (start)').to_dict()['KEGG.organism']
    else:
        loci_in_KEGG = None

    # map to KEGG
    gene_table['KEGG.reaction'] = pd.Series(dtype='str')
    gene_table['KEGG.notes'] = 'updated EC with KEGG: False'
    gene_table = gene_table.progress_apply(lambda row: map_to_KEGG_row(row,loci_in_KEGG), axis=1)

    # save information
    # with KEGG.reaction for adding to model
    # .....................................................................
    # @TODO
    #    is this filter too hard? only EC is a bit risky
    #    how can one get more hits ?
    #    -> best idea would probably be get more EC numbers to check
    # .....................................................................
    print(F'\t\tgenes mapped to KEGG + EC: {len(gene_table.loc[gene_table["KEGG.reaction"].notna()])}')
    gene_table.loc[gene_table['KEGG.reaction'].notna()].to_csv(Path(working_dir,'genes_mapped.csv'), index=False, header=True)
    # without for possible manual curation
    print(F'\t\tgenes NOT mapped to KEGG + EC: {len(gene_table.loc[gene_table["KEGG.reaction"].isna()])}, saved for manual curation')
    gene_table.loc[gene_table['KEGG.reaction'].isna()].to_csv(Path(manual_dir,'genes_no_reaction.csv'), index=False, header=True)

    return Path(working_dir,'genes_mapped.csv')


# mapping to BiGG
# ---------------

# @TODO
def map_BiGG_reactions_row(row:pd.Series, namespace:pd.DataFrame) -> pd.Series:
    """Map a single entry from the table in map_BiGG_reactions() to the BiGG reaction namespace.

    @TODO
        NOTE: only works for cases, where KEGG.reaction in row contains EXACTLY one entry
              in the rare case that multiple reactions belong to one enzyme, they are omitted
              in this search
    
    Args:
        - row (pd.Series): A single row of the table.
        - namespace (pd.DataFrame): The BiGG reaction namespace table.

    Returns:
        pd.Series:  The edited row.
    """
    
    # match by EC number AND KEGG id
    matches = namespace.loc[namespace['EC Number'].str.contains(row['EC number']) & namespace['KEGG Reaction'].str.contains(row['KEGG.reaction'])]

    # case 1 : no matches
    if matches.empty:
        return row

    # case 2 : exactly one match
    elif len(matches) == 1:
        row['bigg_id'] = matches['id'].values[0]

    # case 3 : multiple matches
    #          often due to reaction being in different compartments
    else:
        row['bigg_id'] = ' '.join(matches['id'].values)

    return row


# @TEST : fitted to refinegems
# @CHECK : connections, e.g. input is now a param short 
def map_BiGG_reactions(table_file:str) -> pd.DataFrame:
    """Map the output of map_to_KEGG() to a BiGG namespace file (rewritten-type, see auxilliaries).

    Args:
        - table_file (str): The path to the saved table from running map_to_KEGG().

    Returns:
        pd.DataFrame: The table with an additional column for the mapping to BiGG reactions.
    """

    r_namespace = load_a_table_from_database('bigg_reactions', False)

    table = pd.read_csv(table_file)
    table['bigg_id'] = pd.Series(dtype='str')

    table = table.apply(lambda row: map_BiGG_reactions_row(row,r_namespace), axis=1)

    return table

# ------------
# extent model
# ------------

# @TODO
def isreaction_complete(reac:cobra.Reaction, 
                        exclude_dna:bool=True, exclude_rna:bool=True) -> bool:
    """Check, if a reaction is complete and ready to be added to the model.
    Additionally, it is possible to check for DNA and RNA reations
    and set them to be excluded or included.

    Args:
        - reac (cobra.Reaction): The reaction to be checked.
        - exclude_dna (bool, optional): Tag to include or exclude DNA reactions. 
            Defaults to True.
        - exclude_rna (bool, optional): Tag to include or exclude RNA reactions. 
            Defaults to True.

    Returns:
        bool: True if the check is successful, else false.
    """

    # ................
    # @TODO
    # extendable
    # ................

    # check reaction
    if exclude_dna and 'DNA' in reac.name:
        return False
    if exclude_rna and 'RNA' in reac.name:
        return False

    # check metabolites
    for m in reac.metabolites:
        if m.id == '' or pd.isnull(m.id):
            return False
        if m.name == '' or pd.isnull(m.name):
            return False
        if m.formula == '' or pd.isnull(m.formula):
            return False

    return True


# build genes
# -----------

def add_gene(model:cobra.Model, reaction:str, row:pd.Series, 
             first:bool=False) -> cobra.Model:
    """Add a new gene to a genome-scale metabolic cobra model.

    Args:
        - model (cobra.Model): The model.
        - reaction (str): The reaction id to add the gene to.
        - row (pd.Series): A single row of the output table of map_BiGG_reactions().
        - first (bool, optional):  Shows, if gene is the first gene to be added to the reaction. 
            True if gene is first to be added.
            Defaults to False.

    Returns:
        cobra.Model: Updated model.
    """

    # add gene
    if first or model.reactions.get_by_id(reaction).gene_reaction_rule == '':
        model.reactions.get_by_id(reaction).gene_reaction_rule = row[0]
    else:
        model.reactions.get_by_id(reaction).gene_reaction_rule = model.reactions.get_by_id(reaction).gene_reaction_rule + ' or ' + row[0]

    # add name
    model.genes.get_by_id(row[0]).name = row[1]

    # add annotations
    if not pd.isnull(row[4]):
        model.genes.get_by_id(row[0]).annotation['ncbigene'] = row[4]
    model.genes.get_by_id(row[0]).annotation['ncbiprotein'] = row[2].split('.')[0]
    # note: annotations like sbo, kegg.genes and uniprot missing

    return model


# build metabolites
# -----------------
# @TODO
# @DOCS wrong
# UNDER CONSTRUCTION
def build_metabolite_mnx(metabolite: str, model:cobra.Model, 
                         mnx_chem_prop:pd.DataFrame, mnx_chem_xref:pd.DataFrame, 
                         bigg_metabolites:pd.DataFrame, namespace:pd.DataFrame) -> cobra.Metabolite:
    """Create or retrieve (from model) a metabolite based on its MetaNetX ID.

    Args:
        - metabolite (str): The MetaNetX ID of the metabolite.
        - model (cobra.Model):  The underlying genome-scale metabolic model.
        - mnx_chem_prop (pd.DataFrame): The chem_xref table from MetaNetX.
        - mnx_chem_xref (pd.DataFrame): The chem_prop table from MetaNetX.
        - bigg_metabolites (pd.DataFrame): The BiGG compound namespace table.
        - namespace (pd.DataFrame): The namespace of the model.

    Returns:
        cobra.Metabolite: The retrieved or newly build metabolite.
    """

    metabolite_prop = mnx_chem_prop[mnx_chem_prop['ID']==metabolite]
    metabolite_anno = mnx_chem_xref[mnx_chem_xref['ID']==metabolite]
    model_mnx = [x.annotation['metanetx.chemical'] for x in model.metabolites if 'metanetx.chemical' in x.annotation]

    # fast check if compound already in model
    # ------------------------------------------
    # @TODO ..........................................
        #   currently no checking for compartments
        #   first match will be taken (most often cytosol one)
        #   regardless of the compartment
        #.............................................
    # step 1: check if MetaNetX ID in model
    if metabolite in model_mnx:
        matches = [x.id for x in model.metabolites if 'metanetx.chemical' in x.annotation and x.annotation['metanetx.chemical']==metabolite]

    # step 2: if yes, retrieve metabolite from model
    #  case 1: multiple matches found
        if len(matches) > 1:
            # ................
            # @TODO see above
            # ................
            match = model.metabolites.get_by_id(matches[0])
        #  case 2: only one match found
        else:
            match = model.metabolites.get_by_id(matches[0])

        # step 3: add metabolite
        return match

    # if not, create new metabolite
    # -----------------------------
    else:

        # step 1: create a random metabolite ID
        # ...........................
        # @TODO : compartment problem 
        # - does it have to be in the name?
        # ...........................
        new_metabolite = cobra.Metabolite(create_random_id(model, 'meta','SPECIMEN')) 


        # step 2: add features
        # --------------------
        # @TODO ..........................................
        #   currently no checking for compartments
        #   defaults to c
        #   makes it difficult to add exchange reactions
        #.................................................
        new_metabolite.formula = metabolite_prop['formula'].iloc[0]
        new_metabolite.name = metabolite_prop['name'].iloc[0]
        new_metabolite.charge = metabolite_prop['charge'].iloc[0]
        new_metabolite.compartment = 'c'

        # step 3: add notes
        # -----------------
        new_metabolite.notes['added via'] = 'metanetx.chemical'

        # step 4: add annotations
        # -----------------------
        new_metabolite.annotation['metanetx.chemical'] = metabolite_prop['ID'].iloc[0]
        new_metabolite.annotation['chebi'] = metabolite_prop['reference'].iloc[0].upper()
        if not pd.isnull(metabolite_prop['InChIKey'].iloc[0]):
            new_metabolite.annotation['inchikey'] = metabolite_prop['InChIKey'].iloc[0].split('=')[1]
        for db in ['kegg.compound','metacyc.compound','seed.compound','bigg.metabolite']:
            db_matches = metabolite_anno[metabolite_anno['source'].str.contains(db)]
            if len(db_matches) == 1:
                 new_metabolite.annotation[db] = db_matches['source'].iloc[0].split(':',1)[1]
            elif len(db_matches) > 1:
                new_metabolite.annotation[db] = [m.split(':',1)[1] for m in db_matches['source'].tolist()]

        # if no BiGG was found in MetaNetX, try reverse search in BiGG
        if metabolite in bigg_metabolites['MetaNetX (MNX) Chemical']:
            new_metabolite.annotation['bigg.metabolite'] =  bigg_metabolites[bigg_metabolites['MetaNetX (MNX) Chemical']==metabolite].iloc[0]
        
        # add additional information from bigg if possible    
        if 'bigg.metabolite' in new_metabolite.annotation.keys():
            bigg_information = bigg_metabolites[bigg_metabolites['bigg_id'].str.contains('|'.join(new_metabolite.annotation['bigg.metabolite']))]
            db_id_bigg = {'BioCyc':'biocyc', 'MetaNetX (MNX) Chemical':'metanetx.chemical','SEED Compound':'seed.compound','InChI Key':'inchikey'}
            for db in db_id_bigg:
                info = bigg_information[db].dropna().to_list()
                if len(info) > 0:
                    info = ','.join(info)
                    info = [x.strip() for x in info.split(',')] # make sure all entries are a separate list object
                    new_metabolite.annotation[db_id_bigg[db]] = info

        # step 5: change ID according to namespace
        # ----------------------------------------
        match_id_to_namespace(new_metabolite,namespace)
       
        # step 6: re-check existence of ID in model
        # -----------------------------------------
        # @TODO : check complete annotations? 
        #        - or let those be covered by the duplicate check later on?
        if new_metabolite.id in [_.id for _ in model.metabolites]:
            return model.metabolites.get_by_id(new_metabolite.id)
        
    return new_metabolite

# @TODO
# @DOCS
# UNDER CONSTRUCTION
def build_metabolite_kegg(kegg_id:str, model:cobra.Model, 
                          model_kegg_ids:list[str], bigg_metabolites:pd.DataFrame, 
                          namespace:Literal['BiGG']='BiGG') -> cobra.Metabolite:
    """Create or retrieve (from model) a metabolite based on its KEGG ID.

    Args:
        - kegg_id (str): The KEGG.compound ID of the metabolite in question.
        - model (cobra.Model): The model.
        - model_kegg_ids (list[str]): List of all annotated KEGG Compound IDs in the model.
        - bigg_metabolites (pd.DataFrame):  The BiGG compound namespace table.
        - namespace (Literal[&#39;BiGG&#39;], optional): Namespace of the model. Defaults to 'BiGG'.

    Returns:
        cobra.Metabolite:  The retrieved or newly build metabolite.
    """

    # retrieve KEGG entry for compound
    # --------------------------------
    try:
        kegg_handle = REST.kegg_get(kegg_id)
        kegg_record = [r for r in Compound.parse(kegg_handle)][0]
    except urllib.error.HTTPError:
        print(F'HTTPError: {kegg_id}')
        return cobra.Metabolite()
    except ConnectionResetError:
        print(F'ConnectionResetError: {kegg_id}')
        return cobra.Metabolite()
    except urllib.error.URLError:
        print(F'URLError: {kegg_id}')
        return cobra.Metabolite()

    # ---------------------------------------
    # fast check if compound already in model
    # ---------------------------------------
    # @TODO ..........................................
        #   currently no checking for compartments
        #   first match will be taken (most often cytosol one)
        #   regardless of the compartment
        #.............................................
    # step 1: check via KEGG ID
    if kegg_id in model_kegg_ids:
        matches = [x.id for x in model.metabolites if ('kegg.compound' in x.annotation and x.annotation['kegg.compound'] == kegg_id)]

        # step 2: model id --> metabolite object
        #  case 1: multiple matches found
        if len(matches) > 1:
            match = model.metabolites.get_by_id(matches[0])
        #  case 2: only one match found
        else:
            match = model.metabolites.get_by_id(matches[0])

        # step 3: add metabolite
        return match

    # -----------------------------
    # if not, create new metabolite
    # -----------------------------
    # ...............
    # @TODO
    #     compartment
    # ...............
    else:
        # step 1: create a random metabolite ID
        # -------------------------------------
        # ...........................
        # @TODO : compartment problem 
        # - does it have to be in the name?
        # ...........................
        new_metabolite = cobra.Metabolite(create_random_id(model, 'meta','SPECIMEN')) 

        # step 2: add features
        # --------------------
        # @TODO ..........................................
        #   currently no checking for compartments
        #.............................................
        # set name from KEGG and additionally use it as ID if there is none yet
        if isinstance(kegg_record.name, list):
            if len(kegg_record.name) > 1:
                new_metabolite.name = kegg_record.name[1]
            else:
                new_metabolite.name = kegg_record.name[0]
        else:
            new_metabolite.name = kegg_record.name
        # set compartment
        new_metabolite.compartment = 'c'
        # set formula
        new_metabolite.formula = kegg_record.formula

        # step 3: add notes
        # -----------------
        new_metabolite.notes['added via'] = 'KEGG.compound'

        # step 4: add annotations
        # -----------------------
        new_metabolite.annotation['kegg.compound'] = kegg_id
        db_idtf = {'CAS':'cas','PubChem':'pubchem.compound','ChEBI':'chebi'}
        for db,ids in kegg_record.dblinks:
            if db in db_idtf:
                if len(ids) > 1:
                    new_metabolite.annotation[db_idtf[db]] = ids
                else:
                    new_metabolite.annotation[db_idtf[db]] = ids[0]

        # add additional information from BiGG
        if kegg_id in bigg_metabolites['KEGG Compound']:

            bigg_information = bigg_metabolites[bigg_metabolites['KEGG Compound']==kegg_id]
            if len(bigg_information) > 0:

                new_metabolite.annotation['bigg.metabolite'] = bigg_information['bigg_id'].to_list()

                db_id_bigg = {'BioCyc':'biocyc', 'MetaNetX (MNX) Chemical':'metanetx.chemical','SEED Compound':'seed.compound','InChI Key':'inchikey'}
                for db in db_id_bigg:
                    info = bigg_information[db].dropna().to_list()
                    if len(info) > 0:
                        info = ','.join(info)
                        info = [x.strip() for x in info.split(',')] # make sure all entries are a separate list object
                        new_metabolite.annotation[db_id_bigg[db]] = info

        # step 5: change ID according to namespace
        # ----------------------------------------
        match_id_to_namespace(new_metabolite,namespace)
       
        # step 6: re-check existence of ID in model
        # -----------------------------------------
        # @TODO : check complete annotations? 
        #        - or let those be covered by the duplicate check later on?
        if new_metabolite.id in [_.id for _ in model.metabolites]:
            return model.metabolites.get_by_id(new_metabolite.id)

        return new_metabolite

# @TODO
# @DOCS
# UNDER CONSTRUCTION
def get_metabolites_mnx(model:cobra.Model,equation:str,
                        mnx_chem_xref:pd.DataFrame,mnx_chem_prop:pd.DataFrame,
                        bigg_metabolites:pd.DataFrame, namespace:Literal['BiGG']='BiGG') -> dict[cobra.Metabolite,int]:
    """Based on a given MetaNetX equation and a model, get or
    create metabolite entires in/for the model.

    Args:
        - model (cobra.Model): The GEM.
        - equation (str): The equation from MetaNetX.
        - mnx_chem_xref (pd.DataFrame):  The chem_xref table from MetaNetX.
        - mnx_chem_prop (pd.DataFrame):  The chem_prop table from MetaNetX.
        - bigg_metabolites (pd.DataFrame): The BiGG compound namespace table.
        - namespace (Literal['BiGG'], optional): Namespace of the model. 
            Defaults to 'BiGG'.

    Returns:
        dict[cobra.Metabolite,int]: Dictonary of metabolites and stoichiometric factors.
    """  

    # @TODO ...................................
    #   currently no checking for compartments
    #..........................................

    model_metabolites = [m.formula for m in model.metabolites]
    metabolites = {}
    produced = -1.0
    factor = 0

    for s in equation.split(' '):
        # switch from reactants to products
        if s == '=':
            produced = 1.0
        # found stoichiometric factor
        elif s.isnumeric():
            factor = float(s)
        # skip
        elif s == '+':
            continue
        # found metabolite
        else:
            # get information from MetaNetX
            metabolite, compartment = s.split('@')
            # build or identify metabolite
            new_metabolite = build_metabolite_mnx(metabolite, model, mnx_chem_prop, mnx_chem_xref,bigg_metabolites, namespace)
            # add metabolite
            if new_metabolite.id in [_.id for _ in metabolites]:
                # ......................................................
                # @TODO: 
                #   check if metabolite if both reactant and product
                #   suggests exchange reaction 
                #   -> maybe a good place to change compartment for one?
                #   -> what about name and directions???
                # ......................................................
                try:
                    test = model.metabolites.get_by_id(new_metabolite.id)
                    new_metabolite = new_metabolite.copy()
                    new_metabolite.id = new_metabolite.id + '_i'
                except:
                    new_metabolite.id = new_metabolite.id + '_i'

            metabolites[new_metabolite] = factor * produced

    return metabolites


# @TODO
# @DOCS
# UNDER CONSTRUCTION
def get_metabolites_kegg(model:cobra.Model,equation:str,
                         chem_xref:pd.DataFrame,chem_prop:pd.DataFrame,
                         bigg_metabolites:pd.DataFrame, 
                         namespace:Literal['BiGG']='BiGG') -> dict[cobra.Metabolite,int]:
    """Based on a given KEGG equation and a model, get or
    create metabolite entires in/for the model.

    Args:
        - model (cobra.Model): A GEM.
        - equation (str): The equation from KEGG.
        - chem_xref (pd.DataFrame): The chem_xref table from MetaNetX.
        - chem_prop (pd.DataFrame): The chem_prop table from MetaNetX.
        - bigg_metabolites (pd.DataFrame): The BiGG compound namespace table.
        - namespace (Literal['BiGG'], optional): The namespace of the model. 
            Defaults to 'BiGG'.

    Returns:
        dict[cobra.Metabolite,int]: Dictonary of metabolites and stoichiometric factors.
    """

    # @TODO ...................................
    #   currently no checking for compartments
    #..........................................

    model_metabolites = [m.formula for m in model.metabolites]
    model_kegg_ids = [m.annotation['kegg.compound'] for m in model.metabolites if 'kegg.compound' in m.annotation]
    metabolites = {}
    produced = -1.0
    factor = 1
    mnx_id = ''

    for s in equation.split(' '):
        # switch from reactants to products
        if '=' in s:
            produced = 1.0
        # found stoichiometric factor
        elif s.isnumeric():
            factor = float(s)
        # skip
        elif s == '+':
            continue
        # found metabolite
        else:
            # check if s is a valid ID
            if '(' in s:
                s = s.split('(')[0]
                # ..................................
                # @TODO
                #     known case: DNA(n) --> DNA(n+1)
                #     currently note in brackets gets ignored
                # ..................................
            elif not s.isalnum():
                print('Problem: unknown character in ID inside get_metabolites_kegg() detected.\nPlease contact dev about your problem.')
                sys.exit(1)

            mnx_id = chem_xref[chem_xref['source'] == F'kegg.compound:{s}']
            # case 1:
            #     add metabolite via MetaNetX
            #     -> make sure, only 1 ID match is found (match is unambiguous)
            if len(mnx_id) == 1:
                mnx_id = mnx_id['ID'].item()
                metabolite = build_metabolite_mnx(mnx_id, model, chem_prop, chem_xref, bigg_metabolites, namespace)
            # case 2:
            #     add metabolite via KEGG
            else:
                metabolite = build_metabolite_kegg(s, model, model_kegg_ids, bigg_metabolites, namespace)

            # add new metabolite
            # @TODO : place to check for exchanges?
            if metabolite.id in [_.id for _ in metabolites]:
                try:
                    test = model.metabolites.get_by_id(metabolite.id)
                    metabolite = metabolite.copy()
                    metabolite.id = metabolite.id + '_i'
                except:
                    metabolite.id = metabolite.id + '_i'

            metabolites[metabolite] = factor * produced

    return metabolites


# build reaction
# --------------
#@TODO
# UNDER CONSTRUCTION
def add_reaction(model:cobra.Model,row:pd.Series,
                 reac_xref:pd.DataFrame,reac_prop:pd.DataFrame,
                 chem_xref:pd.DataFrame,chem_prop:pd.DataFrame,
                 bigg_metabolites:pd.DataFrame, 
                 namespace:str='BiGG', 
                 exclude_dna:bool=True, exclude_rna:bool=True) -> cobra.Model:
    """Helper function for extent_model. For a given row of a gene to be added, add the 
    corresponding reactions and metabolites via MetaNetX/KEGG as needed and available.

    Args:
        - model (cobra.Model): The model.
        - row (pd.Series): A row from the table in function extend_model.
        - reac_xref (pd.DataFrame): The MetaNetX reac_xref table.
        - reac_prop (pd.DataFrame): The MetaNetX reac_prop table.
        - chem_xref (pd.DataFrame): The MetaNetX chem_xref table.
        - chem_prop (pd.DataFrame): The MetaNetX chem_prop table.
        - bigg_metabolites (pd.DataFrame): The BiGG compound namespace table.
        - namespace (str, optional): Namespace of the model. 
            Defaults to 'BiGG'.
        - exclude_dna (bool, optional): Flag to exclude the addition of 
            directly DNA-related reactions. 
            Defaults to True.
        - exclude_rna (bool, optional): Flag to exclude the addition of 
            directly RNA-related reactions. 
            Defaults to True.

    Returns:
        cobra.Model: The extended model.
    """

    # create reaction object
    reac = cobra.Reaction(create_random_id(model,'reac','SPECIMEN'))

    # ----------------------------
    # curate reaction via MetaNetX
    # ----------------------------
    # try kegg.reaction --> metanetx.reaction
    if F'kegg.reaction:{row["KEGG.reaction"]}' in list(reac_xref['source']):
        
        # get MetaNetX ID
        met_reac_kegg = reac_xref[reac_xref['source']==F'kegg.reaction:{row["KEGG.reaction"]}']
        met_reac = reac_prop[reac_prop['ID']==met_reac_kegg['ID'].iloc[0]]

        # make sure exactly one entry is parsed
        # @TODO : parallel parsing
        if len(met_reac) > 1:
            print(F'Warning: multiple matches for kegg.reaction {row["KEGG.reaction"]} found. Only first one will be used.')
            met_reac = met_reac.head(1)

        # add name
        # --------
        #     from MetaNetX KEGG description
        reac.name = met_reac_kegg['description'].iloc[0].split('|')[0]

        # add notes
        # ---------
        reac.notes['creation'] = 'via MetaNetX'
        reac.notes['KEGG.information'] = row['KEGG.notes']

        # add metabolites
        # ----------------
        reac.add_metabolites(get_metabolites_mnx(model,met_reac['mnx equation'].iloc[0],chem_xref,chem_prop,bigg_metabolites, namespace))
        #@TODO .............
        #   direction of reaction
        #   ---> current solution:
        #        use one direction only
        # ..................

        # add annotations
        # ---------------
        reac.annotation['ec-code'] = row['EC number']
        reac.annotation['kegg.reaction'] = row['KEGG.reaction']
        reac.annotation['metanetx.reaction'] = met_reac_kegg['ID'].iloc[0]
        met_reac_anno = reac_xref[reac_xref['ID']==met_reac_kegg['ID'].iloc[0]]
        for db in ['metacyc.reaction','seed.reaction','rhea','bigg.reaction']:
            db_matches = met_reac_anno[met_reac_anno['source'].str.contains(db)]
            if len(db_matches) == 1:
                reac.annotation[db] = db_matches['source'].iloc[0].split(':',1)[1]
            elif len(db_matches) > 1:
                reac.annotation[db] = [r.split(':',1)[1] for r in db_matches['source'].tolist()]
            else:
                continue
    
    # if not possible, use information from KEGG only
    # ------------------------
    # curate reaction via KEGG
    # ------------------------
    else:
        
        # retrieve reaction information from KEGG
        reac_kegg = kegg_reaction_parser(row['KEGG.reaction'])

        # add name
        # --------
        #     from KEGG name
        reac.name = reac_kegg['name']

        # add notes
        # ---------
        reac.notes['creation'] = 'via KEGG'
        reac.notes['KEGG.information'] = row['KEGG.notes']

        # add metabolites
        # ----------------
        reac.add_metabolites(get_metabolites_kegg(model,reac_kegg['equation'],chem_xref,chem_prop,bigg_metabolites, namespace))
            #@TODO .............
            #   direction of reaction
            #   ---> current solution:
            #        use one direction only
            # ..................

        # add annotations
        # ---------------
        reac.annotation['ec-code'] = row['EC number']
        reac.annotation['kegg.reaction'] = row['KEGG.reaction']
        for db, identifiers in reac_kegg['db'].items():
            if len(identifiers) == 1:
                reac.annotation[db] = identifiers[0]
            else:
                reac.annotation[db] = identifiers


    # --------------------------------------
    # re-set ID to fit namespace if possible
    # --------------------------------------
    match_id_to_namespace(reac, namespace)

    # ---------------------
    # add reaction to model
    # ---------------------
    
    # if the ID change results in an ID already in the model, use that reaction
    if reac.id in [_.id for _ in model.reactions]:
        print(f'{reac.id} already in model, not added a second time.')
    else:
        # check if reaction is complete
        # and fullfills the requirements / parameters
        if isreaction_complete(reac, exclude_dna, exclude_rna):
            model.add_reactions([reac])
        else:
            print(F'reaction {reac.name} for gene {row["locus_tag"]} could not be completely reconstructed, not added to model.')
            return model

    # --------
    # add gene
    # --------
    # check if gene is already in model
    if row['locus_tag'] in model.genes:
        # if - for whatever reason - gene already in gpr, skip
        if row['locus_tag'] in model.reactions.get_by_id(reac.id).gene_reaction_rule:
            return model
        # create new gpr, if nonexistent
        elif not model.reactions.get_by_id(reac.id).gene_reaction_rule or len(model.reactions.get_by_id(reac.id).gene_reaction_rule) == 0:
            model.reactions.get_by_id(reac.id).gene_reaction_rule = row['locus_tag']
        # add gene to existing gpr
        else:
            model.reactions.get_by_id(reac.id).gene_reaction_rule = model.reactions.get_by_id(reac.id).gene_reaction_rule + ' or ' + row['locus_tag']
    else:
        # add (to) gene reaction rule and curate new gene object
        if not model.reactions.get_by_id(reac.id).gene_reaction_rule or len(model.reactions.get_by_id(reac.id).gene_reaction_rule) == 0:
            model = add_gene(model, reac.id, row, first=True)
        else:
            model = add_gene(model, reac.id, row, first=False)

    return model


# extend model
# ------------
# notes
# @CHECK : connections, e.g. input is now a param short 
def extent_model(table:pd.DataFrame, model:cobra.Model,
                 chem_prop_file:str,chem_xref_file:str,
                 reac_prop_file:str,reac_xref_file:str, 
                 namespace:str='BiGG', 
                 exclude_dna:bool=True, exclude_rna:bool=True) -> cobra.Model:
    """Add reactions, metabolites and genes to a model based on the output of map_to_bigg().

    Args:
        - table (pd.DataFrame): The table with the information to be added to the model.
            Output of :py:func:`map_to_bigg`.
        - model (cobra.Model): The model to extend.
        - chem_prop_file (str): Path to the MetaNetX chem_prop file.
        - chem_xref_file (str): Path to the MetaNetX chem_xref file.
        - reac_prop_file (str): Path to the MetaNetX reac_prop file.
        - reac_xref_file (str): Path to the MetaNetX reac_xref file.
        - namespace (str, optional): Namespace of the model. 
            Defaults to 'BiGG'.
        - exclude_dna (bool, optional): Tag to include or exclude DNA reactions. 
            Defaults to True.
        - exclude_rna (bool, optional): Tag to include or exclude RNA reactions. 
            Defaults to True.

    Returns:
        cobra.Model: The extended model.
    """

    # load MetaNetX database / namespace
    chem_prop = pd.read_csv(chem_prop_file, sep='\t', comment='#', names=['ID','name','reference','formula','charge','mass','InChI','InChIKey','SMILES'])
    chem_xref = pd.read_csv(chem_xref_file, sep='\t', comment='#', names=['source','ID','description'])

    reac_prop = pd.read_csv(reac_prop_file, sep='\t', comment='#', names=['ID','mnx equation','reference','classifs','is_balanced','is_transport'])
    reac_xref = pd.read_csv(reac_xref_file, sep='\t', comment='#', names=['source','ID','description'])

    # load bigg metabolite namespace
    bigg_metabolites = load_a_table_from_database('bigg_metabolites', False)
    bigg_metabolites.rename(columns={'id':'bigg_id'}, inplace=True)
    bigg_metabolites = bigg_metabolites[['bigg_id','universal_bigg_id','name','CHEBI','BioCyc','KEGG Compound','MetaNetX (MNX) Chemical','SEED Compound','InChI Key']]

    # add genes one by one to model
    print('\tAdding genes and if needed reactions and metabolites to model:')
    for row_idx in tqdm(table.index):

        # generate Name -> KEGG.reaction dictionary
        react_dict = get_reaction_annotation_dict(model,'KEGG')
        # generate Name -> BiGG.reaction dictionary
        react_dict_2 = get_reaction_annotation_dict(model,'BiGG')

        # get row in pandas format
        row = table.iloc[row_idx]


        # case 1: BiGG name already in model = reaction in model
        if not pd.isnull(row['bigg_id']) and any((True for _ in row['bigg_id'].split(' ') if _ in react_dict_2.values())):
            # get matching reaction id(s)
            reac_found = [_ for _ in row['bigg_id'].split(' ') if _ in react_dict_2.values()]
            # add genes to all reactions
            for r in reac_found:
                model = add_gene(model, r, row)

        # case 1: KEGG reaction ID in model = reaction probably in model as well
        elif row['KEGG.reaction'] in react_dict.values():
            # get corresponding reaction
            react_found = [_ for _ in react_dict.keys() if row['KEGG.reaction'] == react_dict[_]]
            # add gene to all reactions found
            for r in react_found:
                model = add_gene(model,r,row)

        # case 3: reaction not in model
        #         -> add reaction(s), gene and metabolites if needed
        else:
            # case 3.1: one reaction
            react = row['KEGG.reaction'].split(' ')
            if len(react) == 1:
                model = add_reaction(model,row,reac_xref,reac_prop,chem_xref,chem_prop,bigg_metabolites, namespace, exclude_dna, exclude_rna)

            # case 3.2: multiple reactions
            #           add each reaction separatly, with the same gene for th gene reaction rule
            # note: zero reactions not possible due to previous filtering
            else:

                for r in react:
                    temp_row = row.copy(deep=True)
                    temp_row['KEGG.reaction'] = r
                    model = add_reaction(model,temp_row,reac_xref,reac_prop,chem_xref,chem_prop,bigg_metabolites, namespace, exclude_dna, exclude_rna)

    return model

# main
# ----

# @TEST : fitted to refinegems
# @CHECK : connections, e.g. input is now two params short 
# run the main as function
def run(draft:str, gene_list:str, fasta:str, 
        db:str, dir:str, 
        mnx_chem_prop:str, mnx_chem_xref:str, 
        mnx_reac_prop:str, mnx_reac_xref:str, 
        ncbi_map:str, ncbi_dat:str, 
        id:str='locus_tag', 
        sensitivity:Literal['sensitive', 'more-sensitive', 'very-sensitive','ultra-sensitive']='more-sensitive', 
        coverage:float=95.0, pid:float=90.0, 
        threads:int=2, 
        exclude_dna:bool=True, exclude_rna:bool=True, 
        memote:bool=False):
    """Extend a draft model.

    Using a multitude of information from MetaNetX, KEGG, NCBI and more,
    extend a draft model by yet unadded genes ("missing genes") by mapping them against
    a database using DIAMOND.

    Args:
        - draft (str): Path to the draft model.
        - gene_list (str): Path to a csv file containing information on all the genes found in the annotated genome.
        - fasta (str): Path to the (protein) FASTA file containing the CDS sequences.
        - db (str): Path to the database used for running DIAMOND.
        - dir (str):  Path to the directory for the output (directories).

        - mnx_chem_prop (str): Path to the MetaNetX chem_prop namespace file.
        - mnx_chem_xref (str): Path to the MetaNetX chem_xref namespace file.
        - mnx_reac_prop (str): Path to the MetaNetX reac_prop namespace file.
        - mnx_reac_xref (str): Path to the MetaNetX reac_xref namespace file.

        - ncbi_map (str): Path to the ncbi information mapping file. Optional, but recommended.
        - ncbi_dat (str): Path to the ncbi database information file. Optional, but recommended.

        - id (str, optional): Name of the column of the csv file that contains the entries that were 
            used as gene identifiers in the draft model. 
            Defaults to 'locus_tag'.
        - sensitivity (Literal['sensitive', 'more-sensitive', 'very-sensitive','ultra-sensitive'], optional): Sensitivity mode for DIAMOND blastp run.
            Can be sensitive, more-sensitive, very-sensitive or ultra-sensitive.
            Defaults to 'more-sensitive'.
        - coverage (float, optional): Threshold value for the query coverage for DIAMOND. 
            Defaults to 95.0.
        - pid (float, optional): PID (percentage identity value) to filter the blast hist by. 
            Only hits equal or above the given value are kept.
            Defaults to 90.0.
        - threads (int, optional): Number of threads to be used for DIAMOND. 
            Defaults to 2.

        - exclude_dna (bool, optional): Exclude reactions with DNA in their name when added. 
            Defaults to True.
        - exclude_rna (bool, optional): Exclude reactions with RNA in their name when added. 
            Defaults to True.

        - memote (bool, optional): Use memote on the extended model. 
            Defaults to False.

    Raises:
        ValueError: Unknown sensitivity mode.
    """
    

    if not sensitivity in ['sensitive','more-sensitive','very-sensitive','ultra-sensitive']:
        raise ValueError(F'Unknown sensitivity mode {sensitivity}. Choose from: sensitive, more-sensitive, very-sensitive, ultra-sensitive')
        sys.exit(1)

    print('\nrefinement step 1: extension\n################################################################################\n')


    # -----------------------
    # create output directory
    # -----------------------

    try:
        Path(dir,"step1-extension").mkdir(parents=True, exist_ok=False)
        print(F'Creating new directory {str(Path(dir,"step1-extension"))}')
    except FileExistsError:
        print('Given directory already has required structure.')

    try:
        Path(dir,"manual_curation").mkdir(parents=True, exist_ok=False)
        print(F'Creating new directory {str(Path(dir,"manual_curation"))}')
    except FileExistsError:
        print('Given directory already has required structure.')

    # ----------------------
    # identify missing genes
    # ----------------------

    print('\n# ----------------------\n# identify missing genes\n# ----------------------')

    start = time.time()

    # read in the draft model
    draft = cobra.io.read_sbml_model(Path(draft))
    # identify missing genes and create new FASTA
    missing_genes = identify_missing_genes(gene_list, draft, id, fasta, Path(dir,'step1-extension'))
    print(F'\tmissing genes: {len(missing_genes)}')

    end = time.time()
    print(F'\ttime: {end - start}s')

    # ------------------------------------------------------
    # use DIAMOND to identify homologs for the missing genes
    # ------------------------------------------------------

    print('\n# -----------\n# run DIAMOND\n# -----------')

    start = time.time()

    outname_diamond = Path(dir,'step1-extension','DIAMOND_results.tsv')
    bl = "\\ "
    print(F'\tRunning the following command:')
    print(F'diamond blastp -d {db.replace(" ",bl)} -q {Path(dir.replace(" ",bl),"step1-extension","missing_genes.fasta")} --{sensitivity} --query-cover {coverage} -p {int(threads)} -o {outname_diamond} --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore')
    subprocess.run([F'diamond blastp -d {db.replace(" ",bl)} -q {Path(dir.replace(" ",bl),"step1-extension","missing_genes.fasta")} --{sensitivity} --query-cover {coverage} -p {int(threads)} -o {outname_diamond} --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore'], shell=True)

    end = time.time()
    print(F'\ttime: {end - start}s')

    # -----------------------------------------
    # find best hits for each gene, if possible
    # -----------------------------------------

    print('\n# --------------------\n# filter out best hits\n# --------------------')

    start = time.time()

    # best hit for each query / missing gene, if one exists
    best_hits = find_best_diamond_hits(outname_diamond, pid)
    print(F'\tgenes after DIAMOND with coverage > {coverage}, PID {pid}: {len(best_hits)}')

    end = time.time()
    print(F'\ttime: {end - start}s')

    # --------------------------------------------
    # get further information on the missing genes
    # --------------------------------------------

    print('\n# --------------------------\n# get additional information\n# --------------------------')

    start = time.time()

    # from NCBI
    print('\tRetrieve information from NCBI via accession version number')
    genes_to_add = get_ncbi_info(best_hits, ncbi_map)

    # map to KEGG (gene or ec --> enzyme and reaction)
    print('\tMap to KEGG')
    genes_to_add = map_to_KEGG(genes_to_add, Path(dir,"step1-extension"), Path(dir,"manual_curation"), ncbi_dat)

    # map to BiGG
    print('\tmap information to BiGG namespace via EC number AND KEGG.reaction ID')
    genes_to_add = map_BiGG_reactions(genes_to_add)

    end = time.time()
    print(F'\ttime: {end - start}s')

    # ------------------
    # add genes to model
    # ------------------

    print('\n# ------------------\n# add genes to model\n# ------------------')

    start = time.time()
    # save starting values
    m_before = len(draft.metabolites)
    r_before = len(draft.reactions)
    g_before = len(draft.genes)

    # extent the model
    draft = extent_model(genes_to_add,draft,mnx_chem_prop,mnx_chem_xref,mnx_reac_prop,mnx_reac_xref, exclude_dna, exclude_rna)
    # save it
    name = F'{draft.id}_extended'
    cobra.io.write_sbml_model(draft, Path(dir,'step1-extension',name+'.xml'))

    # find out the differences to before
    r_after = len(draft.reactions)
    print(F'\tnumber of added reactions: {r_after - r_before}')
    m_after = len(draft.metabolites)
    print(F'\tnumber of added metabolites: {m_after - m_before}')
    g_after = len(draft.genes)
    print(F'\tnumber of added genes: {g_after - g_before}')

    end = time.time()
    print(F'\ttime: {end - start}s')

    # ---------------------------------
    # assess model quality using memote
    # ---------------------------------

    if memote:
        memote_path = Path(dir,'step1-extension',name+'.html')
        run_memote(draft, 'html', return_res=False, save_res=memote_path, verbose=True)
