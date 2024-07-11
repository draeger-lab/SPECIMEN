``HQTB`` Configuration File
===========================

Below, the configuration file with the underlying defaults, is shown.

.. code-block:: yaml 
    
    # Configuration file for the SPECIMEN HQTB pipeline

    # Meaning of the default parameters:
    #    The value __USER__ indicates parameters required to be specified by the user
    #    The value USER indicates parameters required only in specific cases

    # Meta info:
    #    model:     USER
    #    organism:  USER
    #    date:      USER
    #    author:    USER

    # Input for the pipeline
    # ----------------------

    # Information about the genome to be used to generate the new model
    subject:
        annotated_genome: __USER__
        full_sequence: __USER__

    # Information about the template model/genome
    template:
        annotated_genome: __USER__
        model: __USER__
        namespace: BiGG

    # Information about the output
    out:
        dir: ./specimen_run/
        name: specimen_model
        memote: False

    # Data(bases) required to run the program
    data:
        # If this parameter is set, assumes that the directory structure from setup
        # is used and uses this path to a directory as the parent folder for the
        # following paths (assumes all data paths are relative ones)
        data_direc: null
        # Required
        diamond: __USER__
        # Needed but potentially downloaded
        mnx_chem_prop: MetaNetX/chem_prop.tsv
        mnx_chem_xref: MetaNetX/chem_xref.tsv
        mnx_reac_prop: MetaNetX/reac_prop.tsv
        mnx_reac_xref: MetaNetX/reac_xref.tsv
        # Optional, but good and manual
        ncbi_map: null
        ncbi_dat: null
        # Optional for directionality control
        biocyc: null
        # Optional:
        #   The pan-core model is used for analysis and if no universal model
        #   is given, also for gapfilling.
        #   If the pan-core model is too small for useful gapfilling, use an
        #   additional universal model for gapfilling.
        #   If none is given gapfilling (and core-pan analysis) is skipped
        universal: null
        pan-core: null

    # Paramters for the single steps of the pipeline
    parameters:
        bidirectional_blast:
            # Default should suffice except special cases
            template_name: null
            input_name: null
            temp_header: null
            in_header: null
            # Can be set by user if wanted, but not necessary
            sensitivity: more-sensitive

        generate_draft_model:
            edit_names: no
            pid: 80.0
            medium: default

        refinement_extension:
            # Default (usually) fine
            id: locus_tag
            # Default fine
            sensitivity: more-sensitive
            # Default alright but good to edit for trying different options
            coverage: 95.0
            pid: 90.0
            # Default almost needed, except for special cases
            exclude_dna: True
            exclude_rna: True

        refinement_cleanup:
            # Default as standard
            check_dupl_reac: True
            check_dupl_meta: default
            remove_unused_meta: False
            remove_dupl_reac: True
            remove_dupl_meta: True
            # Current default means no gapfilling
            media_gap: null

        refinement_annotation:
            # For KEGG pathway annotation
            viaEC: False
            viaRC: False

        refinement_smoothing:
            # Useful
            mcc: skip
            # ECG correction
            egc: null
            # Depend on organism (current: Klebsiella )
            dna_weight_frac: 0.023
            ion_weight_frac: 0.05

        # Validation:
            # Default should suffice

        analysis:
            # Default is currently only option
            pc_based_on: id
            # Can be default but useful to edit
            media_analysis: __USER__ # Edit to fit a default media config file
            test_aa_auxotrophies: True
            # Perform pathway analysis with KEGG
            pathway: True

    # Options for performance
    performance:
        threads: 2
        # For the gapfilling, if iterations and chunk_size are set (not null)
        # use a heuristic for faster performance:
        #     Instead of using all reactions that can be added at once,
        #     run x interations of gapfilling with n-size randomised chunks of reactions
        gapfilling:
            iterations: 3
            chunk_size: 2000
