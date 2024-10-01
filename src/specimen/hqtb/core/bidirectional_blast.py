"""Perform a bidirectional blastp using DIAMOND on an input and a template (annotated genomes).
"""

__author__ = 'Carolin Brune'

################################################################################
# requirements
################################################################################

import pandas as pd
pd.options.mode.chained_assignment = None # editing on copy and saving is elsewhere

import subprocess
import time
import os.path

from Bio import SeqIO
from pathlib import Path
from typing import Literal

################################################################################
# functions
################################################################################

def extract_cds(file: str, name: str, dir: str, collect_info: list, identifier: str) -> str:
    """Extract the CDS from a genbank file (annotated genome). 
    Produces a FASTA-file.

    Args:
        - file (str): 
            File to extract CDS from.
        - name (str): 
            Name of the genome.
        - dir (str): 
            Directory for the ouput.
        - collect_info (list): 
            Feature identifiers to collect information from.
        - identifier (str): 
            Feature identifier to use of the header of the FASTA.

    Returns:
        str: 
            Name of the FASTA-file 
    """
    
    extension = os.path.splitext(os.path.basename(file))[1]


    match extension:

        case '.gbff':
            fasta_name = Path(dir,'FASTA',name+'_prot.fa')
            # check, if CDS were already extracted
            if os.path.isfile(Path(dir,'FASTA',name+'_nucs.fa')) and os.path.isfile(Path(dir,'FASTA',name+'_prot.fa')):
                print(F'CDS for {name} were already extracted')
            else:
                # run extraction process
                info_list = []
                with (
                    open(file) as f,
                    open(Path(dir,'FASTA',name+'_nucs.fa'), "w") as nucs_f,
                    open(Path(dir,'FASTA',name+'_prot.fa'), "w") as prot_f,
                ):

                    for record in SeqIO.parse(f,"genbank"):
                        x = 1
                        if record.features:
                            for feature in record.features:

                                # -----------------
                                # extract sequences
                                # -----------------

                                if feature.type == 'CDS':
                                    # extract nucleotide sequence
                                    nuc_seq = feature.location.extract(record).seq
                                    # extract protein sequence or translate
                                    if 'translation' in feature.qualifiers.keys():
                                        prot_seq = feature.qualifiers['translation'][0]
                                    else:
                                        prot_seq = str(nuc_seq.translate())
                                    # create header
                                    if identifier in feature.qualifiers.keys():
                                        header = F">{feature.qualifiers[identifier][0]}"
                                    else:
                                        header = F">unidentified_{x}"
                                        x += 1
                                    # write to file
                                    nucs_f.write(F'{header}\n{nuc_seq}\n')
                                    prot_f.write(F'{header}\n{prot_seq}\n')

                                    # ---------------------------------
                                    # extract extra feature information
                                    # ---------------------------------
                                    temp_dict = {collect_info[0]:feature.qualifiers[collect_info[0]][0]}
                                    for key in collect_info[1:]:
                                        if key in feature.qualifiers.keys():
                                            temp_dict[key] = feature.qualifiers[key][0]
                                        else:
                                            temp_dict[key] = None
                                    info_list.append(temp_dict)

                info_df = pd.DataFrame.from_dict(info_list)
                print(F'{len(info_df)} CDS were extracted from {name}')
                info_df.to_csv(Path(dir,name+'_info.csv'), index=False)
            return fasta_name

        case '.faa':
            # no CDS extraction needed, only the info-file
            print('faa extension detected. Assuming no CDS extraction is needed.')
            info_df = pd.DataFrame([seq.id for seq in SeqIO.parse(file, 'fasta')], columns=['locus_tag'])
            info_df.to_csv(Path(dir,name+'_info.csv'), index=False)
            fasta_name = file
            return fasta_name

        case _:
            mes = F'Extract_cds: Unkown file extension {extension} for file {file}.'
            raise ValueError(mes)


def create_diamond_db(dir: str, name: str, path: str, threads: int):
    """Create a DIAMOND database for a given protein FASTA file.

    Args:
        - dir (str): 
            Path to the data directory.
        - name (str): 
            Name of the genome/database.
        - path (str): 
            Path to the FASTA-file.
        - threads (int): 
            Number of threads to use.
    """

    # check if database already exists
    if os.path.isfile(Path(dir,'db',name+'.dmnd')):
        print(F'database for {name} already exists')
    else:
        # generate new database using diamond makedb
        print(F'create DIAMOND database for {name} using:')
        print(F'diamond makedb --in {path} -d {str(Path(dir,"db",name+".dmnd"))} -p {int(threads)}')
        start = time.time()
        subprocess.run(["diamond", "makedb", "--in", path, "-d", str(Path(dir,"db",name+".dmnd")), "-p", str(threads)], shell=True)
        end = time.time()
        print(F'\t time: {end - start}s')


def run_diamond_blastp(dir: str, db: str, query: str, fasta_path:str , sensitivity: str, threads: int):
    """Run DIAMOND blastp for a given database name and FASTA - relies on the structure
    created by this file (bidirectional_blast.py).

    Args:
        - dir (str): 
            Parent directory of the place to save the files to.
        - db (str): 
            Name of the genome/database used as the database.
        - query (str): 
            Name of the genome used as the query.
        - fasta_path (str): 
            Path to the FASTA-file containing the CDS.
        - sensitivity (str): 
            Sensitivity mode to use for DIAMOND blastp.
        - threads (int): 
            Number of threads that will be used for running DIAMOND
    """

    outname = Path(dir,'DIAMONDblastp',f'{query}_vs_{db}.tsv')
    # check if file already exists
    if os.path.isfile(outname):
        print(F'file or filename for {query} vs. {db} blastp already exists')
    else:
        # blast file
        print(F'blast {query} against {db} using:')
        print(F'diamond blastp -d {str(Path(dir,"db/",db+".dmnd"))} -q {F"{fasta_path}"} --{sensitivity} -p {int(threads)} -o {outname} --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen ')
        start = time.time()
        subprocess.run(["diamond", "blastp", "-d", str(Path(dir,"db",db+".dmnd")), "-q", fasta_path, "--"+sensitivity, "-p", str(threads), "-o", outname, "--outfmt", "6", "qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qlen"], shell=True)
        end = time.time()
        print(F'\ttime: {end - start}s')


def bdbp_diamond(dir: str, template_name:str, input_name: str, template_path: str, input_path: str, sensitivity='sensitive', threads=2):
    """Perform bidirectional blastp using DIAMOND.

    Args:
        - dir (str): 
            Path to the directory parent to in/out.
        - template_name (str): 
            Name of the template genome.
        - input_name (str): 
            Name of the input genome.
        - template_path (str): 
            Path to the CDS FASTA-file of the template.
        - input_path (str): 
            Path to the CDS FASTA-file of the input.
        - sensitivity (str, optional): 
            Sensitivity mode for DIAMOND. 
            Defaults to 'sensitive'.
        - threads (int, optional): 
            Number of threads to use when running DIAMOND. 
            Defaults to 2.
    """

    # -----------------------
    # create output directory
    # -----------------------
    try:
        Path(dir,"db").mkdir(parents=True, exist_ok=False)
        print(F'Creating new directory {str(Path(dir,"db"))}')
    except FileExistsError:
        print('Given directory already has required structure.')

    try:
        Path(dir,"DIAMONDblastp").mkdir(parents=True, exist_ok=False)
        print(F'Creating new directory {str(Path(dir,"DIAMONDblastp"))}')
    except FileExistsError:
        print('Given directory already has required structure.')

    # ---------------------
    # make DIAMOND database
    # ---------------------

    create_diamond_db(dir, template_name, template_path, threads)
    create_diamond_db(dir, input_name, input_path, threads)

    # ----------------------
    # perform diamond blastp
    # ----------------------

    run_diamond_blastp(dir, template_name, input_name, input_path, sensitivity, threads)
    run_diamond_blastp(dir, input_name, template_name, template_path, sensitivity, threads)


def extract_bestbdbp_hits(tvq: str, qvt: str, name: str, cov:float=0.25):
    """Extract the best directional blastp hits from two tsv files,
    which were generate by :py:func`~specimen.core.bidirectional_blast.bdbp_diamond` 
    generated or similar steps.

    Args:
        - tvq (str): 
            Path to the template vs. query file.
        - qvt (str): 
            Path to the query vs. template file.
        - name (str): 
            Name (path) of the output file.
        - cov (float, optional): 
            Cut-off value for the coverage.
            All hits with coverage < cov will be excluded. 
            Defaults to 0.25.
    """

    # default value for coverage (currently) based on Norsigian et al. 2020

    # read in data and filter for coverage and percentage identity values (PID)
    col_names = ['query_ID', 'subject_ID', 'PID', 'align_len', 'no_mismatch', 'no_gapopen', 'query_start', 'query_end', 'subject_start', 'subject_end','E-value','bitscore','query_length', 'reciprocal']
    # template vs. query
    bh_tvq = pd.read_csv(tvq, sep="\t", names=col_names)
    bh_tvq['coverage'] = bh_tvq['align_len']/bh_tvq['query_length']
    bh_tvq = bh_tvq[bh_tvq['coverage']>=cov]
    # query vs. template
    bh_qvt = pd.read_csv(qvt, sep="\t", names=col_names)
    bh_qvt['coverage'] = bh_qvt['align_len']/bh_qvt['query_length']
    bh_qvt = bh_qvt[bh_qvt['coverage']>=cov]
    # compare dataframes to extract bidirectional best hits
    list_best_hits = list()
    for id in bh_qvt['query_ID'].unique():
        # get best hit from query vs template
        all_hits = bh_qvt[bh_qvt['query_ID'] == id]
        if len(all_hits) == 0:
            continue
        best_hit = all_hits.loc[all_hits['PID'].idxmax()]
        # get best reciprocal hit, if it exists
        all_re_hits = bh_tvq[bh_tvq['query_ID'] == best_hit['subject_ID']]
        if len(all_re_hits) == 0:
            continue
        best_re_hit = all_re_hits.loc[all_re_hits['PID'].idxmax()]
        # add to list of best hits, if they are the same
        if best_re_hit.subject_ID == best_hit.query_ID:
            best_hit['reciprocal'] = 1
            list_best_hits.append(best_hit)
        else:
            best_hit['reciprocal'] = 0
            list_best_hits.append(best_hit)
    # write results into file
    out = pd.DataFrame(list_best_hits, columns=col_names)
    print(F'{len(out)} blast best hits with the template genome')
    print(F'{len(out[out["reciprocal"]==1])} reciprocal and {len(out[out["reciprocal"]==0])} query blast best hits')
    out.to_csv(F"{name}.tsv", index=False, sep="\t")


def run(template:str, input:str, dir:str,
        template_name:str=None, input_name:str=None, 
        temp_header:str='protein_id', in_header:str='locus_tag', 
        threads:int=2, 
        extra_info:list[str]=['locus_tag', 'product', 'protein_id'], 
        sensitivity:Literal['sensitive', 'more-sensitive', 'very-sensitive', 'ultra-sensitive']='more-sensitive'):
    """Run the bidirectional blast on a template and input genome (annotated).

    Args:
        - template (str): 
            Path to the annotated genome file used as a template.
        - input (str): 
            Path to the annotated genome file used as a input.
        - dir (str): 
            Path to the output directory.
        - template_name (str, optional): 
            Name of the annotated genome file used as a template.. 
            Defaults to None.
        - input_name (str, optional): 
            Name of the annotated genome file used as input.. 
            Defaults to None.
        - temp_header (str, optional): 
            Feature qualifier of the gbff (NCBI) / faa (PROKKA) of the template to use as header for the FASTA files.
            If None is given, sets it based on file extension (currently only implemented for gbff and faa).
            Defaults to 'protein_id'.
        - in_header (str, optional): 
            Feature qualifier of the gbff (NCBI) / faa (PROKKA) of the input to use as header for the FASTA files. 
            If None is given, sets it based on file extension (currently only implememted for gbff and faa).
            Defaults to 'locus_tag'.
        - threads (int, optional):  
            Number of threads to be used for DIAMOND. 
            Defaults to 2.
        - extra_info (list[str], optional):  
            List of feature qualifiers to be extracted from the annotated genome files as additional information. 
            Defaults to ['locus_tag', 'product', 'protein_id'].
        - sensitivity (Literal['sensitive', 'more-sensitive', 'very-sensitive', 'ultra-sensitive'], optional): 
            Sensitivity mode for DIAMOND blastp run.. 
            Defaults to 'more-sensitive'.

    Raises:
        - ValueError: Unknown file extension. Please set value for temp_header manually or check file.
        - ValueError: Unknown file extension. Please set value for in_header manually or check file.
        - ValueError: Unknown sensitive mode
    """

    total_time_s = time.time()
    # -----------
    # check input
    # -----------
    if not template_name:
        template_name = os.path.splitext(os.path.basename(template))[0]
    if not input_name:
        input_name = os.path.splitext(os.path.basename(input))[0]

    if not temp_header:
        if '.gbff' == os.path.splitext(template)[1]:
            temp_header = 'protein_id'
        elif '.faa' == os.path.splitext(template)[1]:
            temp_header = 'locus_tag'
        else:
            raise ValueError('Unknown file extension. Please set value for temp_header manually or check file.')

    if not in_header:
        if '.gbff' == os.path.splitext(input)[1]:
            in_header = 'protein_id'
        elif '.faa' == os.path.splitext(input)[1]:
            in_header = 'locus_tag'
        else:
            raise ValueError('Unknown file extension. Please set value for in_header manually or check file.')


    if sensitivity not in ['sensitive','more-sensitive','very-sensitive','ultra-sensitive']:
        raise ValueError(F'Unknown sensitive mode {sensitivity}. Please choose one of the following: sensitive, more-sensitive, very-sensitive, ultra-sensitive')

    # ......................
    # @TODO
    #    when to use protein_id and when to use locus_tag as default (prokka vs ncbi)
    # ......................

    # -------------
    # start program
    # -------------

    print('\nbidirectional blast\n################################################################################\n')

    # ---------------------
    # extract CDS sequences
    # ---------------------

    print('\n# ----------------------\n# extract CDS into FASTA\n# ----------------------\n')

    try:
        Path(dir,"FASTA").mkdir(parents=True, exist_ok=False)
        print(F'Creating new directory {str(Path(dir,"FASTA"))}')
    except FileExistsError:
        print('Given directory already has required structure.')

    print('processing template genome...')
    start = time.time()
    temp_fasta = extract_cds(template, template_name, dir, extra_info, temp_header)
    end = time.time()
    print(F'\ttotal time: {end - start}s')

    print('processing genome of interest...')
    start = time.time()
    inpu_fasta = extract_cds(input, input_name, dir, extra_info, in_header)
    end = time.time()
    print(F'\ttotal time: {end - start}s')


    # ----------------------------
    # perform bidirectional blastp
    # ----------------------------

    print('\n# ----------------------------\n# perform bidirectional blastp\n# ----------------------------\n')

    # reciprical blast - diamond version
    bdbp_diamond(dir, template_name, input_name, temp_fasta, inpu_fasta, sensitivity, threads)

    # ----------------------------
    # find best bidirectional hits
    # ----------------------------

    print('\n# ----------------------------\n# find best bidirectional hits\n# ----------------------------\n')

    start = time.time()
    extract_bestbdbp_hits(Path(dir,'DIAMONDblastp',f'{template_name}_vs_{input_name}.tsv'),Path(dir,'DIAMONDblastp',F'{input_name}_vs_{template_name}.tsv'), Path(dir,F'{input_name}_{template_name}_bbh'))
    end = time.time()
    print(F'\ttotal time: {end - start}s')

    total_time_e = time.time()
    print(F'total runtime: {total_time_e-total_time_s}s')
