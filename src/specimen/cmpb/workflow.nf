nextflow.enable.dsl=2

params.config = "config.yaml" // or your config path

workflow {

    // Step 1: Create draft model with CarveMe if needed
    draft_model = create_draft_model(params.config)

    // Step 2: Validate compartments
    validated_model = validate_compartments(draft_model, params.config)

    // Step 3: CarveMe correction
    carveme_model = carveme_correction(validated_model, params.config)

    // Step 4: Growth test and analysis after draft
    growth_test(carveme_model, params.config, "after_draft")
    analysis(carveme_model, params.config, "after_draft")

    // Step 5: Gapfilling (KEGG, BioCyc, Gene)
    gapfilled_model = gapfill(carveme_model, params.config)

    // Step 6: Growth test and analysis after gapfill
    growth_test(gapfilled_model, params.config, "after_gapfill")
    analysis(gapfilled_model, params.config, "after_gapfill")

    // Step 7: ModelPolisher (optional)
    polished_model = modelpolisher(gapfilled_model, params.config)

    // Step 8: KEGG Pathway Groups
    kegg_annotated_model = kegg_pathway_groups(polished_model, params.config)

    // Step 9: SBOannotator
    sbo_annotated_model = sboannotator(kegg_annotated_model, params.config)

    // Step 10: Duplicates cleanup
    deduped_model = cleanup_duplicates(sbo_annotated_model, params.config)

    // Step 11: Growth test and analysis after duplicate removal
    growth_test(deduped_model, params.config, "after_duplicate_removal")
    analysis(deduped_model, params.config, "after_duplicate_removal")

    // Step 12: Reaction direction check
    direction_checked_model = check_direction(deduped_model, params.config)

    // Step 13: EGCs (energy generating cycles)
    egc_fixed_model = fix_egcs(direction_checked_model, params.config)

    // Step 14: BOF improvement
    bof_model = improve_bof(egc_fixed_model, params.config)

    // Step 15: MCC (CCC)
    mcc_model = run_mcc(bof_model, params.config)

    // Step 16: Save final model
    save_final_model(mcc_model, params.config)

    // Step 17: Final analysis and reports
    final_analysis(mcc_model, params.config)
}

// Example process definitions

process create_draft_model {
    input:
    path config
    output:
    path 'draft_model.xml'
    script:
    """
    python scripts/create_draft_model.py --config $config --output draft_model.xml
    """
}

process validate_compartments {
    input:
    path draft_model
    path config
    output:
    path 'validated_model.xml'
    script:
    """
    python scripts/validate_compartments.py --input $draft_model --config $config --output validated_model.xml
    """
}

// ...repeat for other steps, adapting the script names and arguments...