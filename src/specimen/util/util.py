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

