################################################################################
# requirements
################################################################################

from Bio import SeqIO
import os
import os.path
import pandas as pd
from pathlib import Path
import re
import sys
import subprocess
import time

# further required programs:
#        - DIAMOND, tested with version 0.9.14

################################################################################
# functions
################################################################################

# create a DIAMOND reference database
# -----------------------------------

def create_DIAMOND_db_from_folder(dir, out, name='database', extension='faa', threads=2):
    """Build a DIAMOND database from a folder containing FASTA files.

    :param dir: Path to the directory to search for FASTA files for the database (recursive file search).
    :type dir: string, required
    :param out: Path of the directory of the output.
    :type out: string, required
    :param name: Name of the created database.
        Default is 'database'.
    :type name: string, optional
    :param extension: File extension of the FASTA files (to determine which files to search for).
        Default is 'faa'.
    :type extension: string, optional
    :param threads: Number of threads to use for DIAMOND.
        Default is 2
    :type threads: int, optional
    """

    # check directory ending
    if not dir.endswith("/"):
        dir = dir + "/"
    if not out.endswith("/"):
        out = out + "/"

    # get fasta file names
    fasta_files = Path(dir).rglob(F'*.{extension}')

    # -----------------------
    # combine to 1 FASTA file
    # -----------------------

    # check if folder already has a combined FASTA
    outname_fasta = F'{out}combinded.faa'
    save = True
    if os.path.isfile(outname_fasta):
        print('A combined.faa files already exists in the given folder.')
        answer = None
        while answer not in ('overwrite','use','end'):
            answer = input('Do you want to [overwrite] or [use] the file or [end] the programm?: ')
            if answer == 'use':
                save = False
            elif anwer == 'end':
                sys.exit()
            elif answer == 'overwrite':
                save = True
            else:
                print('Enter [overwrite] or [use] or [end] (without brackets)')


    # save all in new (combined) FASTA
    if save:
        with open(outname_fasta, 'w') as out:
            for f in fasta_files:
                SeqIO.write(SeqIO.parse(f, 'fasta'), out,'fasta')

    # -------------------------
    # generate DIAMOND database
    # -------------------------

    outname_dnmd = F'{out}{name}'
    bl = "\\ "
    print(F'running the following command:\ndiamond makedb --in {outname_fasta.replace(" ",bl)} -d {outname_dnmd.replace(" ",bl)} -p {threads}')
    subprocess.run([F'diamond makedb --in {outname_fasta.replace(" ",bl)} -d {outname_dnmd.replace(" ",bl)} -p {threads}'], shell=True)
    print('finished')


# create a NCBI mapping file for the database
# -------------------------------------------

def get_info_GenBank_Record(file_path):
    """Retrieves a table containg information about the following qualifiers from a
    Genbank file: ['protein_id','locus_tag','db_xref','old_locus_tag','EC_number'].

    :param file_path: Path to the Genbank (.gbff) file.
    :type  file_path: str
    :returns:         A table containing the information above.
    :rtype:           pd.DataFrame, columns= ['ncbi_accession_version', 'locus_tag_ref','old_locus_tag','GeneID','EC number']
    """

    temp_table = pd.DataFrame(columns=['ncbi_accession_version', 'locus_tag_ref','old_locus_tag','GeneID','EC number'])
    attributes = ['protein_id','locus_tag','old_locus_tag','db_xref','EC_number']

    for record in SeqIO.parse(gbff_name,"genbank"):
        if record.features:
            for feature in record.features:
                if feature.type == 'CDS':
                    temp_list = []
                    for a in attributes:
                        if a in feature.qualifiers.keys():
                            temp_list.append(feature.qualifiers[a][0])
                        else:
                            temp_list.append('-')
                    temp_table.loc[len(temp_table)] = temp_list

    # reformat
    pat = re.compile(r'\D')
    temp_table['GeneID'] = [pat.sub('', x) for x in temp_table['GeneID']]

    return temp_table


def create_NCBIinfo_mapping(dir, out, extension='gbff'):
    """Create a NCBI information mapping file from a folder containing e.g. gbff files.

    :param dir: Path to the directory for the recursive file search for the mapping.
    :type dir: string
    :param out: Path of the directory for the output.
    :type out: string
    :param extension: Name of the file extension to be searched.
        Default is gbff, and currently it is advised to leave it at that.
    :type extension: string
    """

    # check directory ending
    if not dir.endswith("/"):
        dir = dir + "/"

    # get gbff file names
    gbff_files = Path(dir).rglob(F'*.{extension}')
    gbff_num = len(list(Path(dir).rglob(F'*.{extension}')))

    # -----------------------------------------------
    # extract information and create the mapping file
    # -----------------------------------------------

    with open(out, 'w') as out_file:

        file_counter = 1

        for gbff_name in gbff_files:

            print(F'Parsing {file_counter}/{gbff_num}: {gbff_name}')
            # retrieve information
            info = get_info_GenBank_Record(gbff_name)
            # save
            info.to_csv(out_file, header=True, index=False)
            # go to next
            file_counter += 1


# rewrite the BiGG namespaces into tables
# ---------------------------------------

def separate_db_links_reaction(row):
    """Separate the database links in the column of the same name
    from a reaction BiGG namespace file row.

    :param row: one row of the table
    :type  row: pandas object, a pd.DataFrame row
    :returns:   The row with new columns for certain database links (EC number, BioCyc, KEGG, MNX, SEED)
    :rtype:     pandas object, a pd.DataFrame row
    """

    values = {'EC Number': [], 'BioCyc': [], 'KEGG Reaction': [], 'MetaNetX (MNX) Equation': [], 'SEED Reaction': []}
    if isinstance(row['database_links'], str):
        links = row['database_links'].split(';')
        for link in links:
            key, value = link.split(':',1)
            key = key.strip()
            value = value.rsplit('/',1)[1].strip()
            if key in values.keys():
                values[key].append(value)

    for k in values.keys():
        row[k] = ', '.join(values[k])

    return row


def rewrite_reactions(file, out):
    """Rewrites or reformates a given BiGG reaction namespace TXT file to the following columns:
    bigg_id, name, reaction_string, EC number, BioCyc, MetaNetX (MNX) Equation, SEED Reaction

    :param file: Path of the input file.
    :type  file: string
    :param out:  Path of the output file.
    :type  out:  string
    """

    # read in the txt
    data = pd.read_csv(file, sep='\t')
    # remove model list and old bigg ids
    data.drop(columns=['model_list', 'old_bigg_ids'], axis=1, inplace=True)
    # add new columns for the database links
    data = data.apply(separate_db_links_reaction, axis=1)
    data.drop(columns=['database_links'], axis=1, inplace=True)
    # save the reformatted data
    data.to_csv(out, sep='\t', index=False, header=True)


def separate_db_links_metabolite(row):
    """Separate the database links in the column of the same name
    from a metabolite BiGG namespace file row.

    :param row: one row of the table
    :type  row: pandas object, a pd.DataFrame row
    :returns:   The row with new columns for certain database links (CHEBI, BioCyc, KEGG, MNX, SEED, InChI)
    :rtype:     pandas object, a pd.DataFrame row
    """

    values = {'CHEBI': [], 'BioCyc': [], 'KEGG Compound': [], 'MetaNetX (MNX) Chemical': [], 'SEED Compound': [], 'InChI Key': []}
    if isinstance(row['database_links'], str):
        links = row['database_links'].split(';')
        for link in links:
            key, value = link.split(':',1)
            key = key.strip()
            value = value.rsplit('/',1)[1].strip()
            if key in values.keys():
                values[key].append(value)

    for k in values.keys():
        row[k] = ', '.join(values[k])

    return row


def rewrite_metabolites(file, out):
    """Rewrites or reformates a given BiGG metabolites namespace TXT file to the following columns:
    bigg_id, name, reaction_string, EC number, BioCyc, MetaNetX (MNX) Equation, SEED Reaction

    :param file: Path of the input file.
    :type  file: string
    :param out:  Path of the output file.
    :type  out:  string
    """

    # read in the txt
    data = pd.read_csv(file, sep='\t')
    # remove model list and old bigg ids
    data.drop(columns=['model_list', 'old_bigg_ids'], axis=1, inplace=True)
    # add new columns for the database links
    data = data.apply(separate_db_links_metabolite, axis=1)
    data.drop(columns=['database_links'], axis=1, inplace=True)
    # save the reformatted data
    data.to_csv(out, sep='\t', index=False, header=True)


def write_BiGG_namespace_to_table(input, out=None, type='reactions'):
    """Rewrite a BiGG namespace into a table format.

    :param input: The input BiGG namespace txt-file (its path). Can be for reactions or metabolites.
    :type input: string
    :param out: Path to the output file. Default is the name of the input file with the additional tag '_rewritte'.
    :type out: string
    :param type: Specifies if the input file is for reactions or metabolites.
        Can either be 'reactions' (default) or 'metabolites'.
    :type type: string
    """

    print('\nrewrite BiGG txt files\n################################################################################\n')

    if out == None:
        out = os.path.dirname(input) + '/' + os.path.splitext(os.path.basename(input))[0] + '_rewritten.tsv'

    # -------------
    # start program
    # -------------

    match type:

        case 'reactions':
            start = time.time()
            print('\trewriting reaction namespace ...')
            rewrite_reactions(input, out)
            end = time.time()
            print(F'\ttime: {end - start}s')

        case 'metabolites':
            start = time.time()
            print('\trewriting metabolite namespace ...')
            rewrite_metabolites(input, out)
            end = time.time()
            print(F'\ttime: {end - start}s')
        case _:
            raise ValueError(F'Unknown option for namespace type: {type}')
