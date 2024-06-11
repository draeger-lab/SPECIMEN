CMPB Configuration File
=======================

Below, the configuration file with the underlying defaults, is listed.

.. code-block:: yaml 
    
    # Configuration file for the SPECIMEN CMPB pipeline

    # Explaination for default parameters:
    #    with the value __USER__ are required to be specified by the user
    #    with the value USER are required only under specific cases

    # meta info:
    #    model:     USER
    #    organism:  USER
    #    date:      USER
    #    author:    USER

    # input for the pipeline
    # ----------------------
    input:
    modelpath: NULL            # optional, path to a model. 
                                # If not given, runs CarveMe
    annotated_genome: __USER__ # required, path to the annotated genome file
    mediapath: __USER__        # path to a media config to tests growth with

    # general options
    # ---------------
    general:
    dir: './'                  # Path/Name of a directory to save output to
    colours: 'YlGn'            # set the colour scheme for the plots
                                # should be a valid matplotlib continuous color palette
    namespace: BiGG            # Namespace to use for the model
                                # Possible identifiers, currently: BiGG
    save_all_models: True      # save all models (models for each step)
    memote_always_on: False    # run memote after every step
    stats_always_on: False     # calculate the model statistics after every step
    # below are options used by multiple steps
    refseq_gff: USER           # Path to RefSeq GFF file: Required for gap analysis with 'KEGG'. 
                                # Can be optionally provided for cm-polish.
    kegg_organism_id: USER     # KEGG ID of the organism: Required for gap analysis with KEGG.
                                # Can be optionally provided for cm-polish.
    
    # part-specific options
    # ---------------------

    # polish a CarveMe model
    #    only neccessary, if the mode will or has been build with CarveMe
    #    will only be used, if model is indeed a CarveMe model
    cm-polish:
    email: USER              # User Mail to use for Entrez 
    protein_fasta: USER      # optional, except for is_lab_strain: True. 
                            # The path to the protein FASTA used to create the CarveMe model.
    is_lab_strain: False     # whether the users strain originates from a lab 
                            # Needs to be set to ensure that protein IDs get the 'bqbiol:isHomologTo' qualifier
                            # & to set the locus_tag to the ones obtained by the annotation

    # gapfilling, optional
    gapfilling:
    ### Automatic gap filling ###
    # All parameters are required for all db_to_compare choices except biocyc_files which is not required for 'KEGG'
    gap_fill_params:
        db_to_compare: USER  # One of the choices KEGG|BioCyc|KEGG+BioCyc 
        biocyc_files: 
        - USER # Path to TXT file containing a SmartTable from BioCyc with the columns 'Accession-2' 'Reaction of gene'
        - USER  # Path to TXT file containing a SmartTable with all reaction relevant information (*)
        - USER  # Path to TXT file containing a SmartTable with all metabolite relevant information (+)
        - USER  # Path to protein FASTA file used as input for CarveMe (Required to get the protein IDs from the locus tags)
    # (*) 'Reaction' 'Reactants of reaction' 'Products of reaction' 'EC-Number' 'KEGG Reaction' 'MetaNetX' 'Reaction-Direction' 'Spontaneous?'
    # (+) 'Compound' 'Object ID' 'Chemical Formula' 'InChI-Key' 'ChEBI'
    gap_fill_file: NULL # Path to Excel file with which gaps in model should be filled
    # Either obtained by running gapfilling/Created by hand with the same structure as the result file from gapfilling
    # Example Excel file to fill in by hand: data/modelName_gapfill_analysis_date_example.xlsx

    # add KEGG pathways as groups, optional
    kegg_pathway_groups: True  # decide, whether to run this or not

    # resolve duplicates
    duplicates:
    # three possible option for the resolvement of duplicates for the following model entities:
    # - check:  check for duplicates and simply report them
    # - remove: check for and remove duplicates from the model (if possible)
    # - skip:   skip the resolvement 
    reactions: remove
    metabolites: remove
    # additional remove unused metabolites (reduces possible knowledge base)
    remove_unused_metabs: False

    # BOFdat / Biomass objective function
    BOF:
    run_bofdat: False
    # if BOFdat should be run, 
    # fill out the params below
    bofdat_params:
        full_genome_sequence: USER  # whole genome sequence
        dna_weight_fraction: USER   # DNA weight fraction for the organism
        weight_fraction: USER       # Ezyme/ion weight fractions for the organism


