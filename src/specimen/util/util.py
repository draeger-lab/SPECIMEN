################################################################################
# requirements
################################################################################


import os
import os.path
import pandas as pd
import re
import subprocess
import sys

from Bio import SeqIO
from pathlib import Path
from typing import Literal

# further required programs:
#        - DIAMOND, tested with version 0.9.14

################################################################################
# variables
################################################################################

################################################################################
# functions
################################################################################

# create a DIAMOND reference database
# -----------------------------------

def create_DIAMOND_db_from_folder(dir:str, out:str, name:str='database',
                                  extension:str='faa', threads:int=2):
    """Build a DIAMOND database from a folder containing FASTA files.

    Args:
        - dir (str): 
            Path to the directory to search for FASTA files for 
            the database (recursive file search).
        - out (str): 
            Path of the directory of the output.
        - name (str, optional): 
            Name of the created database. 
            Defaults to 'database'.
        - extension (str, optional): 
            File extension of the FASTA files 
            (to determine which files to search for). 
            Defaults to 'faa'.
        - threads (int, optional): 
            Number of threads to use for DIAMOND. 
            Defaults to 2.
    """

    # get fasta file names
    fasta_files = Path(dir).rglob(F'*.{extension}')

    # -----------------------
    # combine to 1 FASTA file
    # -----------------------

    # check if folder already has a combined FASTA
    outname_fasta = Path(out,'combinded.faa')
    save = True
    if os.path.isfile(outname_fasta):
        print('A combined.faa files already exists in the given folder.')
        answer = None
        while answer not in ('overwrite','use','end'):
            answer = input('Do you want to [overwrite] or [use] the file or [end] the programm?: ')
            if answer == 'use':
                save = False
            elif answer == 'end':
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

    outname_dnmd = Path(out,name)
    bl = "\\ "
    print(F'running the following command:\ndiamond makedb --in {outname_fasta.replace(" ",bl)} -d {outname_dnmd.replace(" ",bl)} -p {threads}')
    subprocess.run([F'diamond makedb --in {outname_fasta.replace(" ",bl)} -d {outname_dnmd.replace(" ",bl)} -p {threads}'], shell=True)
    print('finished')


# create a NCBI mapping file for the database
# -------------------------------------------

def get_info_GenBank_Record(file_path:str) -> pd.DataFrame:
    """Retrieves a table containg information about the following qualifiers from a
    Genbank file: ['protein_id','locus_tag','db_xref','old_locus_tag','EC_number'].

    Args:
        - file_path (str): 
            Path to the Genbank (.gbff) file.

    Returns:
        pd.DataFrame: 
            A table containing the information above.
            Has the following  columns= ['ncbi_accession_version', 'locus_tag_ref','old_locus_tag','GeneID','EC number'].
    """

    temp_table = pd.DataFrame(columns=['ncbi_accession_version', 'locus_tag_ref','old_locus_tag','GeneID','EC number'])
    attributes = ['protein_id','locus_tag','old_locus_tag','db_xref','EC_number']

    for record in SeqIO.parse(file_path,"genbank"):
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


def create_NCBIinfo_mapping(dir:str, out:str, extension:Literal['gbff']='gbff'):
    """Create a NCBI information mapping file from a folder containing e.g. gbff files.

    Args:
        - dir (str): 
            Path to the directory for the recursive file search for the mapping.
        - out (str): 
            Path of the directory for the output.
        - extension (Literal['gbff'], optional): 
            Name of the file extension to be searched.
            Default is gbff, and currently it is advised to leave it at that.
            Defaults to 'gbff'.
    """

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

