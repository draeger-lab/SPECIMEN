"""Perform a bidirectional blastp using DIAMOND on an input and a template (annotated genomes).
"""

__author__ = 'Carolin Brune'

################################################################################
# requirements
################################################################################

import pandas as pd
pd.options.mode.chained_assignment = None # editing on copy and saving is elsewhere
from Bio import SeqIO
from pathlib import Path
import os.path
import time
import subprocess

################################################################################
# functions
################################################################################

def extract_cds(file: str, name: str, dir: str, collect_info: list, identifier: str):
    """Extract the CDS from a genbank file (annotated genome). Produces a FASTA
    file.

    :param file:          filename to extract CDS from
    :type  file:          string
    :param dir:           directory for the generated files
    :type  dir:           string
    :param name:          name of the genome
    :type  name:          string
    :param collect_info:  feature identifiers to collect information from
    :type  collect_info:  list
    :param identifier:    feature identifier to use as a header for the FASTA files
    :type  identifier:    string
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
            exit(F'Unkown file extension {extension} for file {file}.')


def create_diamond_db(dir: str, name: str, path: str, threads: int):
    """Create a DIAMOND database for a given protein FASTA file.

    :param dir:      path to the data directory
    :type  dir:      string
    :param name:     name of the genome
    :type  name:     string
    :param path:     path to the FASTA file
    :type  path:     string
    :param threads:  number of threads to be used
    :type  threads:  int
    """
    # check if database already exists
    if os.path.isfile(Path(dir,'DIAMONDdb',name+'.dmnd')):
        print(F'database for {name} already exists')
    else:
        # generate new database using diamond makedb
        # @TODO: check if commands run under different OS
        bl = "\\ "
        print(F'create DIAMOND database for {name} using:')
        print(F'diamond makedb --in {path.replace(" ",bl)} -d {str(Path(dir,"db",name)).replace(" ",bl)} -p {int(threads)}')
        start = time.time()
        subprocess.run([F'diamond makedb --in {path.replace(" ",bl)} -d {str(Path(dir,"db",name)).replace(" ",bl)} -p {int(threads)}'], shell=True)
        end = time.time()
        print(F'\t time: {end - start}s')


def run_diamond_blastp(dir: str, db: str, query: str, fasta_path:str , sensitivity: str, threads: int):
    """Run DIAMOND blastp for a given database name and FASTA - relies on the structure
    created by this file (bidirectional_blast.py).

    :param dir:         parent directory of the place to save the files to.
    :type  dir:         string
    :param db:          name of the genome used as the database.
    :type  db:          string
    :param query:       name of the genome used as the query.
    :type  query:       string
    :param fasta_path:  path of the fasta file (for the CDS)
    :type  fasta_path:  string
    :param sensitivity: sensitivity mode that should be used for blastp.
    :type  sensitivity: string
    :param threads:     number of threads should be used for running DIAMOND
    :type  threads:     int
    """

    outname = Path(dir,'DIAMONDblastp',f'{query}_vs_{db}.tsv')
    # check if file already exists
    if os.path.isfile(outname):
        print(F'file or filename for {query} vs. {db} blastp already exists')
    else:
        # blast file
        bl = "\\ "
        print(F'blast {query} against {db} using:')
        print(F'diamond blastp -d {str(Path(dir,"db/",db+".dmnd")).replace(" ",bl)} -q {F"{fasta_path}".replace(" ",bl)} --{sensitivity} -p {int(threads)} -o {outname.replace(" ",bl)} --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen ')
        start = time.time()
        subprocess.run([F'diamond blastp -d {str(Path(dir,"db",db+".dmnd")).replace(" ",bl)} -q {F"{fasta_path}".replace(" ",bl)} --{sensitivity} -p {int(threads)} -o {outname.replace(" ",bl)} --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen '], shell=True)
        end = time.time()
        print(F'\ttime: {end - start}s')


def bdbp_diamond(dir: str, template_name:str, input_name: str, template_path: str, input_path: str, sensitivity='sensitive', threads=2):
    """Perform bidirectional blastp using DIAMOND.

    :param dir:           path to the working directory
    :type  dir:           string
    :param template_name: name of the template genome
    :type  template_name: string
    :param input_name:    name of the input genome
    :type  input_name:    string
    :param template_path: path of the template genome fasta file (for the CDS)
    :type  template_path: string
    :param input_path:    path of the input genome fasta file (for the CDS)
    :type  input_path:    string
    :param sensitivity:   sensitivity mode of DIAMOND, default is mid-sensitive
    :type  sensitivity:   string
    :param threads:       number of threads to be used for computation, default is 2
    :type  threads:       int
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


def extract_bestbdbp_hits(tvq: str, qvt: str, name: str, cov=0.25):
    """Extract the best directional blastp hits from two tsv files,
    which were generate by bdbp_diamond() generated or similar steps.

    :param tvq:  path to the template_vs_query file.
    :type  tvq:  string
    :param qvt:  path to the query_vs_template file.
    :type  qvt:  string
    :param name: name (path) of the output file.
    :type  name: string
    :param cov:  cutoff value for the coverage. The default if 0.25, which means all hits with coverage <0.25 will be excluded.
    :type  cov:  float
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


def run(template, input, dir,template_name=None, input_name=None, temp_header='protein_id', in_header='locus_tag', threads=2, extra_info=['locus_tag', 'product', 'protein_id'], sensitivity='more-sensitive'):
    """Run the bidirectional blast on a template and inpute genome (annotated).

    :param template: Path to the annotated genome file used as a template.
    :type template: string
    :param input: Path to the annotated genome file used as a input.
    :type input: string
    :param dir: Path to the output directory.
    :type dir: string
    :param template_name: Name of the annotated genome file used as a template.
    :type template_name: string, optional
    :param input_name: Name of the annotated genome file used as input.
    :type input_name: string, optional
    :param temp_header: Feature qualifier of the gbff (NCBI) / faa (PROKKA) of the template to use as header for the FASTA files.
        Default is None. If None is given, sets it based on file extension (currently only implememted for gbff and faa).
    :type temp_header: string, optional
    :param in_header: Feature qualifier of the gbff (NCBI) / faa (PROKKA) of the input to use as header for the FASTA files.
        Default is None. If None is given, sets it based on file extension (currently only implememted for gbff and faa).
    :type in_header: string, optional
    :param threads: Number of threads to be used for DIAMOND.
        Default is 2.
    :type threads: int, optional
    :param extra_info: List of feature qualifiers to be extracted from the annotated genome files as additional information.
        Default is ['locus_tag', 'product', 'protein_id'].
    :type extra_info: list, optional
    :param sensitivity: Sensitivity mode for DIAMOND blastp run.
        Can be sensitive, more-sensitive, very-sensitive or ultra-sensitive.
        Default is 'more-sensitive'.
    :type sensitivity: string
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
