``CMPB`` Configuration File
===========================

Below, the configuration file with the underlying defaults, is displayed.

.. code-block:: yaml 
    

    # Configuration file for the SPECIMEN CMPB workflow

    # Meaning of the default parameters:
    #    The value __USER__ indicates parameters required to be specified by the user
    #    The value USER indicates parameters required only in specific cases
    #    To avoid warnings, set parameters you do not use to null or NULL

    # Meta info:
    #    model:     USER
    #    organism:  USER
    #    date:      USER
    #    author:    USER

    # Input for the pipeline
    # ----------------------
    input:
    modelpath: NULL            # Optional, path to a model.
                               # If not given, runs CarveMe 
    mediapath: __USER__        # Path to a media config to test growth with

    # General options
    # ---------------
    general:
    dir: './'                  # Path/Name of a directory to save output to
    colours: 'YlGn'            # Set the colour scheme for the plots
                               # should be a valid matplotlib continuous color palette
    namespace: BiGG            # Namespace to use for the model
                               # Possible identifiers, currently: BiGG
    save_all_models: True      # Save a model per step
    memote_always_on: False    # Run MEMOTE after every step
    stats_always_on: False     # Calculate the model statistics after every step

    # Options used by multiple steps of the workflow
    refseq_gff: USER           # Path to RefSeq GFF file: 
                                # Can be optionally provided for cm-polish.
    kegg_organism_id: USER     # KEGG ID of the organism: Required for gap analysis with 'KEGG'.
                                # Can be optionally provided for cm-polish.
    protein_fasta: USER        # Required, if used for CarveMe or GeneGapFiller. Optional
                                # for cm-polish except for 'is_lab_strain: True'.
                                # The path to the protein FASTA used to create the CarveMe model.
                                # For more information, please refer to the documentation.

    tech-resources:
    email: USER   # User Mail to use for Entrez (accessing NCBI).
    threads: 2    # Number of threads available for tools like DIAMOND.

    # Part-specific options
    # =====================

    # Build a model using CarveMe
    # ---------------------------
    carveme:
        # CarveMe requires protein_fasta under general to be set instead of modelpath
        # if CarveMe should be run, 
        # fill out the params below
        gram: USER      # Choose either grampos or gramneg, depending on the Gram-test
                        # resilts of your organism

    # Polish a CarveMe model
    #    Only neccessary, if the model will or has been build with CarveMe
    #    Will only be used, if model is indeed a CarveMe model
    cm-polish:
    is_lab_strain: False     # Whether the users strain originates from a lab
                            # Needs to be set to ensure that protein IDs get the 'bqbiol:isHomologTo' qualifier
                            # & to set the locus_tag to the ones obtained by the annotation
                            # (Warning: Might cause issues if annotatione was not performed with NCBI PGAP!)

    # Filling gaps, optional
    # ----------------------
    gapfilling:
    
        ########### general options ###########
        # parameters, that apply to all the gap filling algorithmns
        idprefix: 'CMPB'            # prefix to use for fantasy IDs, if IDs for 
                                    # the namespace do not exist.
        formula-check: 'existence'  # When checking, if a metabolite can be added to the model
                                    # also check the formula. For more information about
                                    # available options, please refer to the docs of 
                                    # the function isreaction_comlete().
        exclude-dna: True           # Exclude reactions containing 'DNA' in their name
                                    # from being added to the model.
        exclude-rna: True           # Exclude reactions containing 'RNA' in their name
                                    # from being added to the model.
        
        ########## enable algorithms ##########
        # via KEGG ...................
        # requires KEGG organism ID to be set
        KEGGapFiller: False   # activate gap filling via GFF
        # via BioCyc .................
        BioCycGapFiller: False        # Activate gap filling via BioCyc.
        BioCycGapFiller parameters:   
            gene-table: USER            # Path to a gene smart table file from BioCyc.
            reacs-table: USER           # Path to a reactions smart table from BioCyc.
            gff: USER                   # Path to a GFF file of the genome of the model.
        # via GFF ....................
        GeneGapFiller: False              # Activate gap filling via GFF
        GeneGapFiller parameters: 
            gff: USER                     # Path to a gff file (does not have to be the RefSeq).
                                          # Needs to be from the same genome the model was build on.
            swissprot-dmnd: USER          # Path to the SwissProt DIAMOND database file.
            swissprot-mapping: USER       # Path to the SwissProt mapping file (against EC / BRENDA)
            check-NCBI: False             # Enable checking NCBI accession numbers for EC numbers - time costly.
            sensitivity: 'more-sensitiv'  # Sensitivity option for the DIAMOND run.
            coverage: 90.0                # Coverage (parameter for DIAMOND).
            percentage identity: 90.0     # Percentage identity threshold value for accepting
                                          # matches found by DIAMOND as homologous.


    # Add KEGG pathways as groups, optional
    # -------------------------------------
    kegg_pathway_groups: True

    # Resolve duplicates
    # ------------------
    duplicates:
    # Three possible options for the resolvement of duplicates for the following model entities:
    # - check:  Check for duplicates and simply report them
    # - remove: Check for and remove duplicates from the model (if possible)
    # - skip:   Skip the resolvement
    reactions: remove
    metabolites: remove
    # Additionally, remove unused metabolites (possibly reduces knowledge-base)
    remove_unused_metabs: False

    # Finding and solvong Energy Generating Cycles (EGCs)
    # ---------------------------------------------------
    EGCs:
    solver: NULL # solver gives the algorithm to use for solving EGCs
                # if NULL, only searches for EGCs without trying to solve them
                # options include: greedy

    # BOFdat / Biomass objective function
    # -----------------------------------
    BOF:
    run_bofdat: False
    # if BOFdat should be run, 
    # fill out the params below
    bofdat_params:
        full_genome_sequence: USER  # Whole genome sequence
        dna_weight_fraction: USER   # DNA weight fraction for the organism
        weight_fraction: USER       # Enzyme/ion weight fractions for the organism