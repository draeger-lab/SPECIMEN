"""Entry points to the code from the command line.
"""

__author__ = 'Carolin Brune'

################################################################################
# requirements
################################################################################

import specimen
import click

################################################################################
# entry points
################################################################################

###############################
# unchanged
###############################

@click.group()
@click.help_option("--help", "-h")
@click.version_option()
def cli():
    """SPECIMEN - strain specific metabolic modelling

    This tools provides scripts and functions to build, curate and analyse
    a strain-specific GEM based on a high quality template model and additional
    database files.
    """

# -----
# setup
# -----
@cli.group()
def setup():
    """Setup tools, folder structure and more for running the program.
    """

# get a config file
# -----------------
@setup.command()
@click.option('--filename', '-f', default='config.yaml', type=str, show_default=True, help='Name (path) to the save the config file under.')
@click.option('--type', '-t', default='htqb-basic', type=click.Choice(['htqb-basic', 'htqb-advanced','htqb-defaults','media','cmpb']), help='Type of config file to download. Either a more detailed one for advanced usage or a basic one for beginners with less options.')
def config(filename,type):
    """Download a configuration file (.yaml).

    Download a configuration file to edit for running a complete
    workflow or setting the definitions for a set of media.
    """
    specimen.util.set_up.download_config(filename, type)


# setup data directory structure
# ------------------------------
@setup.command()
@click.argument('pipeline', type=click.Choice(['hqtb','high-quality template based', 
                                               'cmpb', 'carveme modelpolisher based']))
@click.option('--dir','-d', type=str, default='./data/', show_default=True, help='Name/path to the directory create subdirectories in.')
@click.option('--chunk-size', '-s', type=int, default=2048, show_default=True, help=' Size of the chunks of data while downloading.')
def data_structure(pipeline, dir, chunk_size):
    """Create a data directory and download basic databases.

    Creates the 'ideal' data directory structure and directly downloads
    the MetaNetX and BiGG data files.

    PIPELINE is the type of pipeline for which the data structure should be build.
    """
    specimen.util.set_up.build_data_directories(pipeline, dir, chunk_size)

#################
# hqtb pipeline #
#################
@cli.group()
def hqtb():
    """Pipeline for GEM curation based on a high-quality template."""

# run complete pipeline from config
# ---------------------------------
@hqtb.command()
@click.argument('config', type=str)
def run_pipeline(config):
    """Run the complete pipeline based on a config file.

    CONFIG is the path to the configuration file to read the parameters from.
    """
    specimen.workflow.run_complete(config)


# run complete pipeline from config and folder (run multiple  times)
# ------------------------------------------------------------------
@hqtb.command()
@click.argument('config', type=str)
@click.option('-d','--directory', default='', type=str, help='Path to the (parent) directory that contains the folders if the subject input files.')
def run_wrapper(config,directory):
    """Run the complete pipeline multiple times based on a config file
    and a folder. The folder should contain subfolders with the subject files
    (annotated and full genome).

    CONFIG is the path to the configuration file to read the parameters from.
    """
    specimen.workflow.wrapper_pipeline(config, parent_dir=directory)


# run bidirectional blast
# -----------------------
@hqtb.command()
@click.argument('template', type=str)
@click.argument('input', type=str)
@click.option('--template-name', type=str, help='Name of the annotated genome file used as a template.')
@click.option('--input-name', type=str, help='Name of the annotated genome file used as a input.')
@click.option('--temp-header', type=str, default='protein_id', help='Feature qualifier of the gbff of the template to use as header for the FASTA files')
@click.option('--in-header', type=str, default='locus_tag', help='Feature qualifier of the gbff of the input to use as header for the FASTA files')
@click.option('--dir', '-d', type=str, default='./01_bidirectional_blast/', help='Path to output directory. Will create a new one, if given path does not exist.')
@click.option('--threads', '-t', type=int, default=2, help='Number of threads to be used.')
@click.option('--sensitivity', '-s', type=click.Choice(['sensitive','more-sensitive','very-sensitive','ultra-sensitive']), default='sensitive', help='Sensitivity mode for DIAMOND blastp run. Can be sensitive, more-sensitive, very-sensitive or ultra-sensitive. Default is sensitive.')

def bdb(template, input, template_name, input_name, temp_header, in_header, dir, threads, sensitivity):
    """Step 1 of the pipeline: Perform bidirectional blast on a TEMPLATE and an INPUT annotated genome.

    TEMPLATE is an annotated genome file (path) that is used for comparison.
    INPUT is an annotated genome file (path) that will be compared to TEMPLATE
    """
    specimen.core.bidirectional_blast.run(template, input, dir,template_name, input_name, temp_header, in_header, threads, extra_info=['locus_tag', 'product', 'protein_id'], sensitivity=sensitivity)

###############################
# below untested, but corrected
###############################

# generafte draft
# ---------------
@hqtb.command()
@click.argument('template', type=str)
@click.argument('bpbbh', type=str)
@click.option('--dir', type=str, default='./02_generate_draft_model', help='Path to the output directory.')
@click.option('--edit-names', type=click.Choice(['no', 'dot-to-underscore']), default='no', show_default=True, help='Edit the identifier of the FASTA files to match the names in the model.')
@click.option('--pid', type=float, default=80.0, show_default=True, help='Threshold value (percentage identity) for determining, if a gene is counted as present or absent')
@click.option('--name', type=str, default=None, help='Name of the output model, will be taken from the file name if not specified.')
@click.option('--medium', type=str, default='default', help='Set the medium for the new model. If not set, will use the one from the template. If given the keyword "exchanges", will open all exchange reaction and use them as the medium.')
@click.option('--namespace','--nsp',type=click.Choice(['BiGG']),default='BiGG',help='Namespace of the model.')
@click.option('--memote', is_flag=True, default=False, help='Run Memote on the generated draft model.')
def draft(template, bpbbh, dir, edit_names, pid, name, medium, nsp, memote):
    """Step 2 of the pipeline: Generate a draft model from a blastp best hits tsv file and a template model.

    TEMPLATE is the path (string) to the template model.\n
    BPBBH is the path (string) to the  BLASTp bidirectional best hits (step 1).
    """
    specimen.core.generate_draft_model.run(template, bpbbh, dir, edit_names, 
                                           pid, name, medium, nsp, memote)


# refinement
# ----------
@hqtb.group()
def refinement():
    """Step 3 of the workflow: Refinement of the model
    """

@refinement.command()
@click.option('--draft', type=str, required=True, help='Path to the draft model.')
@click.option('--gene-list', '-g', type=str, required=True, help='Path to a csv file containing information on all the genes found in the annotated genome.')
@click.option('--fasta', '-f', type=str, required=True, help="Path to the (protein) FASTA file containing the CDS sequences")
@click.option('--db', '--database', type=str, required=True, help="string, path to the database used for running DIAMOND.")
# MetaNetX links
@click.option('--mnx-chem-prop', type=str, required=True, help='Path to the MetaNetX chem_prop namespace file.')
@click.option('--mnx-chem-xref', type=str, required=True, help='Path to the MetaNetX chem_xref namespace file.')
@click.option('--mnx-reac-prop', type=str, required=True, help='Path to the MetaNetX reac_prop namespace file.')
@click.option('--mnx-reac-xref', type=str, required=True, help='Path to the MetaNetX reac_xref namespace file.')
# NCBI mapping options
@click.option('--ncbi-map', type=str, required=False, help='Path to the ncbi information mapping file. Optional, but recommended.')
@click.option('--ncbi-dat', type=str, required=False, help='Path to the ncbi database information file. Optional, but recommended.')
# DIAMOND and other params
@click.option('--dir','-d', type=str, default='./refinement/', help='Path to the directory for the output (directories)')
@click.option('--id', '-i', type=str, default='locus_tag', help='Name of the column of the csv file that contains the entries that were used as gene identifiers in the draft model.')
@click.option('--sensitivity', '-s', type=click.Choice(['sensitive','more-sensitive','very-sensitive','ultra-sensitive']), default='sensitive', help='Sensitivity mode for DIAMOND blastp run. Default is sensitive.')
@click.option('--coverage', '-c', type=float, default=80.0, help="Threshold value for the query coverage for DIAMOND. Default is 80.0.")
@click.option('--pid', type=float, default=95.0, help='PID (percentage identity value) to filter the blast hist by. Default is 90.0, only hits equal or above the given value are kept.')
@click.option('--threads', '-t', type=int, default=2, help='Number of threads to be used.')
@click.option('--include-dna', is_flag=True, default=False, help='Include reactions with DNA in their name when added (developer information: True == excluded).')
@click.option('--include-rna', is_flag=True, default=False, help='Include reactions with RNA in their name when added (developer information: True == excluded).')
@click.option('--memote', is_flag=True, default=False, help='Use memote on the extended model.')
def extension(draft, gene_list, fasta, db, dir,
    mnx_chem_prop, mnx_chem_xref, mnx_reac_prop, mnx_reac_xref,
    ncbi_map, ncbi_dat,
    id, sensitivity,
    coverage, pid, threads,
    include_dna, include_rna,
    memote):
    """Refinement step 1: Extend the model.

    The following options are required:
    draft, gene_list, fasta, db, dir,
    mnx_chem_prop, mnx_chem_xref, mnx_reac_prop, mnx_reac_xref
    """
    specimen.core.refinement.extension.run(draft, gene_list, fasta, db, dir,
        mnx_chem_prop, mnx_chem_xref, mnx_reac_prop, mnx_reac_xref,
        ncbi_map, ncbi_dat,
        id, sensitivity,
        coverage, pid, threads,
        not include_dna, not include_rna,
        memote)


@refinement.command()
@click.argument('model', type=str)
@click.option('--dir','-d', type=str, default='./refinement/', help='Path to the directory for the output (directories)')
# check directionality
@click.option('--biocyc_db', type=str, required=False, help='Path to the BioCyc (MetaCyc) database information file (for reactions). Optional, but recommended. Necessary for checking directionality')
# check duplicates
@click.option('--check_dupl_reac', '--cdr', is_flag=True, default=False, help='Check for duplicate reactions.')
@click.option('--check_dupl_meta', '--cdm', default='default', type=click.Choice(['default', 'skip', 'exhaustive']), help='Check for duplicate metabolites. Can "default" (starting point MetaNetX), exhaustive (iterate over all annotations as starting points) or "skip".')
@click.option('--objective_function', '--of', type=str, default='Growth', help='Name, ID of the objective function of the model. Default is "Growth".')
@click.option('--remove_dupl_meta', '--rdm', is_flag=True, default=False, help='Option for removing/replacing duplicate metabolites.')
@click.option('--remove_unused_meta', '--rum', is_flag=True, default=False, help='Option for removing unused metabolites from the model. Only used when cdm is not skipped.')
@click.option('--remove_dupl_reac', '--rdr', is_flag=True, default=False, help='Option for removing duplicate reaction from the model.')
# perform gapfilling
@click.option('--universal', '-u', required=False, type=str, help='Path to a universal model containing reactions used for gapfilling.')
@click.option('--media-path', '--mp', required=False, type=str, default=None, help='Path to a media config to use for gapfilling.')
@click.option('--namespace','--nsp', required=False, type=click.Choice(['BiGG']), help='Namespace to use for the model.')
@click.option('--growth_threshold', '-gt', default=0.05, show_default=True, type=float, help='Threshold value for a model to be considered growing.')
@click.option('--iterations', '-i', type=int, default=3, show_default=True, help='Number of iterations for the gapfilling. If 0 is passed, uses full set of reactions instead of heuristic.')
@click.option('--chunk_size', type=int, default=10000, show_default=True, help='Number of reactions to be tested simultaniously if using the heuristic version of gapfilling. If this is 0, heuristic will not be applied.')
# evaluate with memote
@click.option('--memote', is_flag=True, default=False, help='Use memote on the extended model.')
def cleanup(model,
    dir,
    biocyc_db,
    check_dupl_reac,
    check_dupl_meta,
    objective_function,
    remove_unused_meta,
    remove_dupl_reac,
    remove_dupl_meta,
    universal,
    mp,
    nsp,
    growth_threshold,
    iterations,
    chunk_size,
    memote):
    """Refinement step 2: cleanup

    Includes (all steps are optional):\n
    - directionality check,\n
    - completetion of BioCyc/MetaCyc annotations,\n
    - duplicate removal,\n
    - gapfilling

    MODEL is the path to the model to perform the this refinement step on.
    Ideally in the format of this workflow or the results might differ.
    """
    specimen.core.refinement.cleanup.run(model,
        dir,
        biocyc_db,
        check_dupl_reac,
        check_dupl_meta,
        objective_function,
        remove_unused_meta,
        remove_dupl_reac,
        remove_dupl_meta,
        universal,
        mp,
        nsp,
        growth_threshold,
        iterations,
        chunk_size,
        memote)
    
###############################
# unchanged
###############################

@refinement.command()
@click.argument('model', type=str)
@click.option('--dir','-d', type=str, default='./refinement/', help='Path to the directory for the output (directories)')
@click.option('--kegg-via-ec', '--via-ec','--ec', is_flag=True, default=False, help='Try to map EC numbers to KEGG pathway, if KEGG reaction cannot be mapped directly.')
@click.option('--kegg-via-rc', '--via-rc','--rc', is_flag=True, default=False, help='Try to map KEGG reaction class to KEGG pathway, if KEGG reaction cannot be mapped directly.')
@click.option('--memote', is_flag=True, default=False, help='Use memote on the extended model.')
def annotation(model,dir,kegg_via_ec,kegg_via_rc,memote):
    """Refinement step 3: Annotation

    Further annotate the model.
    Includes:
    - SBO annotation using SBOannotataor
    - adding KEGG pathways based on KEGG reaction ID (+EC,+RC optionally)

    MODEL is the path to the model to be annotated.
    """
    specimen.core.refinement.annotation.run(model,
                                            dir,
                                            kegg_viaEC=kegg_via_ec,
                                            kegg_viaRC=kegg_via_rc,
                                            memote=memote)


###############################
# below untested, but corrected
###############################

@refinement.command()
@click.argument('model',type=str)
@click.option('--genome', '-g', required=True, type=str, help='Path to the genome FASTA (e.g. .fna) file of your genome.')
@click.option('--dir', '-d', default='./refinement/', type=str, help='Path to a directory for the output.')
# additional arguments for EGC solving
@click.option('--egc-solver','--egc', required=False, type=click.Choice([None,'greedy']), default=None, multiple=False, help='String sets the type of solver to use to solve EGCs. Otherwise just reports existing EGCs.')
@click.option('--namespace','--nsp',default='BiGG', type=click.Choice(['BiGG']),help='Namespace of the model.')
# additional arguments for MCC
@click.option('--mcc', default='skip', type=click.Choice(['apply','extra','skip']), help='Option to perform MassChargeCuration on the model. Can be used directly on model or as extra information. Choices are "apply","extra" and "skip". Deafult is "skip".')
# additional arguments for BOF
@click.option('--dna_weight_frac', default=0.023, type=float, help='DNA macromolecular weight fraction for your organism. Default is 0.023 for Klebsiella based on Liao et al.')
@click.option('--ion_weight_frac', default=0.05, type=float, help='weight fraction for the coenzymes and ions. Default is 0.05 based on the default of BOFdat.')
# additional arguments for memote
@click.option('--memote', is_flag=True, default=False, help='Use memote on the extended model.')
def smoothing(model, genome, dir, mcc, dna_weight_frac, ion_weight_frac, egc, namespace, memote):
    """Refinement step 4: Smoothing

    Further refine the model by (optionally) using MCC,
    checking for cycles and using BOFdat.

    MODEL is the path to the model that is to b refined.\n
    Further required is a genome FASTA file of the genome the model was build on.
    """
    specimen.core.refinement.smoothing.run(genome, model, dir, mcc, 
                                           egc,
                                           namespace,
                                           dna_weight_frac, ion_weight_frac, 
                                           memote)


#######################
     # below OK #
#######################

# validation
# ----------
@hqtb.command()
@click.argument('model',type=str)
@click.option('--dir', '-d', default='./validation/', type=str, help='Path to a directory for the output.')
@click.option('--run-test', '-t', multiple=True, default=['all'], help='define, which tests should be run. Current possibilities are "all" and "cobra"')
def validation(model,dir,run_test):
    """Step 4 of the pipeline: Validate the model.

    MODEL is the path to the model to be validated.
    """
    if 'all' in run_test:
        specimen.core.validation.run(dir, model, tests=None, all=True)
    else:
        specimen.core.validation.run(dir, model, tests=run_test, all=False)



# analysis
# --------
@hqtb.command()
# in / out
@click.argument('model', type=str)
@click.option('--dir', default='./analysis/', type=str, help='Path to a directory for the output.')
# pan core analysis
@click.option('--pan-core-comparison', '--pcc', type=str, default='id', help='Option on which feature the comparison of pan-core model and model should should be based on.\nDefault is "id".')
@click.option('--pan-core-model', '--pcm', required=False, type=str, help='Path to a pan-core model.')
# growth analysis
@click.option('--namespace','-n',required=False, type=click.Choice(['BiGG']), multiple=False, default='BiGG', help='Namespace used by the given model.')
@click.option('--media-path','--mp',required=False, type=str, default=None, help='Path to a media config file. Enables growth analysis if given.')
@click.option('--test-aa-auxotrophies', '--taa', is_flag=True, default=False, help='Option to test media/model for auxotrophies.')
@click.option('--pathway', '--pathway-analysis', is_flag=True, default=False, help='Option to perform a pathway analysis using KEGG pathway identifiers.')
def analysis(model,
        dir,
        pcm,
        pcc,
        n,
        mp,
        test_aa_auxotrophies,
        pathway):
    """Step 5 of the pipeline: Analyse the final model.

    Includes a statistical analysis and optional a
    pan-core as well as a growth analysis.

    MODEL is the path to the model to be analysed.
    """
    specimen.core.analysis.run(model_path=model, 
                               dir=dir, 
                               media_path=mp, 
                               namespace=n,
                               pc_model_path=pcm, 
                               pc_based_on=pcc, 
                               test_aa_auxotrophies=test_aa_auxotrophies, 
                               pathway=pathway)


#################
# cmpb pipeline #
#################

@cli.group()
def cmpb():
    """Pipeline for GEM curation based on CarveMe model and ModelPolisher."""

@cmpb.command()
@click.argument('config',type=click.Path(exists=True))
def run(config):
    """Run the pipeline for GEM curation based on a CarveMe model using a config file.

    CONFIG is the path to the config file.
    """
    
    specimen.cmpb.run(config)
