#!/usr/bin/env python
"""Running the CarveMe-ModelPolisher-based (CMPB) workflow."""

__author__ = "Tobias Fehrenbach, Famke Baeuerle, Gwendolyn O. Döbel and Carolin Brune"

# @DISCUSSION modelname, authorinitials, strainid and organism are required if a model is provided. Should we keep this?
################################################################################
# requirements
################################################################################

from datetime import date
import logging
import os
import pandas as pd
from pathlib import Path
import time
from typing import Union
import subprocess
import sys

from cobra import Reaction, Model
from libsbml import Model as libModel

from cobra.io.sbml import _sbml_to_model
from refinegems.analysis import growth
from refinegems.classes.medium import load_media
from refinegems.classes.reports import ModelInfoReport, SBOTermReport
from refinegems.utility.util import test_biomass_presence
from refinegems.curation.biomass import check_normalise_biomass
from refinegems.curation.charges import correct_charges_modelseed
from refinegems.curation.curate import resolve_duplicates, check_direction, polish_model, prune_mass_unbalanced_reacs, add_ECO_terms
from refinegems.curation.miriam import change_all_qualifiers, polish_annotations
from refinegems.classes.gapfill import KEGGapFiller, BioCycGapFiller, GeneGapFiller, KEGGapFiller_via_homolog
from refinegems.curation.pathways import set_kegg_pathways, kegg_pathway_analysis
from refinegems.utility.entities import (
    resolve_compartment_names,
    are_compartment_names_valid,
)
from refinegems.utility.connections import (
    run_memote,
    perform_mcc,
    adjust_BOF,
    run_SBOannotator,
    run_ModelPolisher,
    run_DIAMOND_blastp
)
from refinegems.utility.io import load_model, convert_cobra_to_libsbml, write_model_to_file
from refinegems.developement.decorators import implement
from refinegems.classes import egcs

from ..util.set_up import save_cmpb_user_input, validate_config, build_data_directories

################################################################################
# setup logging
################################################################################
# general logging
logger = logging.getLogger(__name__)

################################################################################
# functions
################################################################################

# @TODO Add debug switches to test growth on every step (Maybe also on full medium?) -> Maybe as entry point?
# @ASK Provide the draft model directly from CarveMe only in debug mode or always? -> Is always provided already.
# dev notes:
#   in the run function: current_model means the cobrapy model,
#   while current_libmodel means the libsbml model
def run(configpath: Union[str, None] = None):
    """Run the CarveMe-ModelPolisher-based (CMPB) workflow.

    Args:
        - configpath (Union[str,None]):
            The Path to a configuration file. If none given, prompts the user to enter the needed
            information step by step.
            Defaults to None.
    """

    def between_save(
        model: Union[Model, libModel],
        model_dir: Union[None, str] = None,
        label: Union[None, str] = None,
        only_modelpath: Union[None, str] = None
    ): 
        """Helper function to save the model in between the different steps of the workflow

        Args:
            - model (Union[Model,libModel]):
                The current model, either loaded with libsbml or COBRApy.
            - model_dir (Union[None,str], optional):
                Output model directory. Defaults to None.
            - label (Union[None,str], optional):
                Name of the step. Defaults to None.
            - only_modelpath (Union[None,str], optional):
                Path to save the model to, if only one model will be saved.
                If set, ignores directory and label.
                Defaults to None.
        """
        name = model.id if isinstance(model, Model) else model.getId()
        if not only_modelpath:
            write_model_to_file(
                model, str(Path(model_dir, f"{name}_{label}.xml"))
            )
        else:
            write_model_to_file(model, str(only_modelpath))

    def between_growth_test(model: Model, cfg: dict, step: str):
        """Helper function for :py:func:`~specimen.cmpb.workflow.run`.
        Run a growth test on a model with the options set in the config file.

        Args:
            - model (Model):
                The cobra.Model for testing the growth.
            - cfg (dict):
                The loaded configuration file.
            - step (str):
                Descriptive string for the current step of the pipeline (important for saving the output).
        """
 
        # try to set objective to growth
        growth_func_list = test_biomass_presence(model)
        if growth_func_list:
            # independently of how many growth functions are found, the first one will be used
            model.objective = growth_func_list[0]
            # simulate growth on different media
            growth_report = growth.growth_analysis(
                model,
                cfg["input"]["mediapath"],
                namespace=cfg["general"]["namespace"],
                retrieve="report",
            )

            # Generate path for report
            growth_dir = Path(cfg["general"]["dir"], "cmpb_out", "misc", "growth", step)
            try:
                Path(growth_dir).mkdir(parents=True, exist_ok=False)
            except FileExistsError:
                logger.warning(
                    f"Given directory {growth_dir} already exists. High possibility of files being overwritten."
                )
            
            growth_report.save(
                growth_dir,
                color_palette=config["general"]["colours"],
            )
        else:
            mes = f"No growth/biomass function detected, growth simulation for step {step} will be skipped."
            logger.warning(mes)

    def between_analysis(model: Model, cfg: dict, step: str):
        """Helper function for :py:func:`~specimen.cmpb.workflow.run`.
        Run a memote and/or stats test on a model with the options set in the config file.

        Args:
            - model (Model):
                The cobra.Model for testing.
            - cfg (dict):
                The loaded configuration file.
            - step (str):
                Descriptive string for the current step of the pipeline (important for saving the output).
        """
        # optional analysis
        if cfg["general"]["memote_always_on"]:
            run_memote(
                model,
                "html",
                return_res=False,
                save_res=Path(
                    cfg["general"]["dir"], "cmpb_out", "misc", "memote", f"{step}.html"
                ),
                verbose=False,
            )
        if cfg["general"]["stats_always_on"]:
            report = ModelInfoReport(model)

            # Generate path for report
            stats_dir = Path(cfg["general"]["dir"], "cmpb_out", "misc", "stats", step)
            try:
                Path(stats_dir).mkdir(parents=True, exist_ok=False)
            except FileExistsError:
                logger.warning(
                    f"Given directory {stats_dir} already exists. High possibility of files being overwritten."
                )
                
            report.save(
                Path(stats_dir),
                color_palette=config["general"]["colours"],
            )

    # setup phase
    #############

    # load config
    # -----------
    if not configpath:
        config = save_cmpb_user_input()
    else:
        config = validate_config(configpath, "cmpb")

    # set up model name
    # -----------------
    if (
        config["general"]["authorinitials"] is not None
        and config["general"]["organism"] is not None
        and config["general"]["strainid"] is not None
    ):
        modelname = (
            "i"
            + config["general"]["organism"]
            + str(config["general"]["strainid"])
            + config["general"]["authorinitials"]
            + str(date.today().year).removeprefix("20")
        )
    elif config["general"]["modelname"] is not None:
        modelname = config["general"]["modelname"]
    else:
        logger.info(
            "No values given for the standard name for a model. Default name will be used."
        )
        modelname = "model_" + str(date.today().year).removeprefix("20")

    # Set up directory structure names
    PARENT_DIR = config["general"]["dir"]
    CMPB_OUT_DIR = Path(PARENT_DIR, 'cmpb_out')
    MODEL_DIR = Path(CMPB_OUT_DIR, 'models')
    MISC_DIR = Path(CMPB_OUT_DIR, 'misc')

    # Get cmodel path if only one model should be saved
    only_modelpath = None
    if not config["general"]["save_all_models"]:
        only_modelpath = Path(CMPB_OUT_DIR, f"{modelname}.xml")

    # create directory structure
    # --------------------------
    build_data_directories("cmpb", PARENT_DIR)

    # create log
    # ----------
    today = date.today().strftime("%Y%m%d")
    log_file = Path(CMPB_OUT_DIR, "logs", f"specimen_cmpb_{str(today)}.log")
    log_file.unlink(missing_ok=True) # Remove file in case it already exists
    handler = logging.handlers.RotatingFileHandler(
        log_file, mode="w", backupCount=10, encoding="utf-8", delay=0
    )
    handler.setFormatter(
        logging.Formatter(
            "{levelname} \t {name} \t {message}",
            style="{",
        )
    )
    logger.addHandler(handler)

    # redirect cobrapy logging
    cobralogger = logging.getLogger("cobra")
    cobralogger.addHandler(handler)
    cobralogger.propagate = False

    # redirect refineGEMs logging
    rglogger = logging.getLogger("refinegems")
    rglogger.addHandler(handler)
    rglogger.propagate = False

    total_time_s = time.time()
    
    # CarveMe
    #########
    if not config["input"]["modelpath"]:
        
        logger.info("Creating a draft model with CarveMe ...")
        step_start = time.time()
        
        out_modelpath = only_modelpath if only_modelpath else Path(MODEL_DIR, modelname+'.xml')
        if config["carveme"]["refseq"] is not None:
            os.system(
                f"carve --refseq {config['carveme']['refseq']} --solver scip -o {out_modelpath}"
            )
        elif (
            config["carveme"]["gram"] == "grampos"
            or config["carveme"]["gram"] == "gramneg"
        ):
            os.system(
                f"carve {config['general']['protein_fasta']} --solver scip -u {config['carveme']['gram']} -o {out_modelpath}"
            )
        else:
            os.system(
                f"carve {config['general']['protein_fasta']} --solver scip -o {out_modelpath}"
            )
        config["input"]["modelpath"] = out_modelpath
        
        step_end = time.time()
        logger.info(f"CarveMe\truntime: {step_end-step_start}s\n")

    # Set model path for pipeline
    current_modelpath = config["input"]["modelpath"]

    # load & validate compartments of the loaded model
    ##################################################
    
    logger.info("Validating compartments ...")
    step_start = time.time()
    
    # load into libsbml
    current_libmodel = load_model(str(current_modelpath), "libsbml")

    # reload into cobrapy & validate compartments
    current_model = _sbml_to_model(current_libmodel.getSBMLDocument())
    # check, if compartment names are valid
    if not are_compartment_names_valid(current_model):
        # if not, attempt to fix them
        resolve_compartment_names(current_model)
        # save validated model
        between_save(
            current_model, MODEL_DIR, "validated_comps", only_modelpath
        )

        # Transform model into libmodel for CarveMe correction
        current_libmodel = convert_cobra_to_libsbml(current_model, 'notes')
        
    step_end = time.time()
    logger.info(f"Compartment validation\truntime: {step_end-step_start}s\n")

    # CarveMe correction
    ####################

    # check, if input is a CarveMe model
    if ("CarveMe" in current_libmodel.getNotesString()) or ("CarveMe" in current_libmodel.getAnnotationString()):
        
        logger.info("Performing CarveMe correction ...")
        step_start = time.time()
        
        # Set up directory for report output
        current_sub_dir = Path(MISC_DIR, "model_correction")
        try:
            Path(current_sub_dir).mkdir(parents=True, exist_ok=False)
            print(f"Creating new directory {current_sub_dir}")
        except FileExistsError:
            logger.warning(
                f"Given directory {current_sub_dir} already exists. High possibility of files being overwritten."
            )

        # Polish model
        current_libmodel = polish_model(
            current_libmodel,
            id_db=config["general"]["namespace"],
            gff_paths=[config["general"]["gff"]],
            email=config["tech-resources"]["email"],
            lab_strain=config["cm-polish"]["is_lab_strain"],
            kegg_organism_id=config["general"]["kegg_organism_id"],
            outpath=str(current_sub_dir),
        )
        # rg correct charges
        current_libmodel, mult_charges_dict = correct_charges_modelseed(
            current_libmodel
        )
        if mult_charges_dict:
            mult_charges_tab = pd.DataFrame.from_dict(mult_charges_dict, orient="index")
            mult_charges_tab.to_csv(
                Path(MISC_DIR, "reac_with_mult_charges.tsv"), sep="\t"
            )

        # save model
        between_save(
            current_libmodel, MODEL_DIR, "after_CarveMe_correction", only_modelpath
        )
        
        step_end = time.time()
        logger.info(f"CarveMe correction\truntime: {step_end-step_start}s\n")

    # growth test
    # -----------
    current_model = _sbml_to_model(current_libmodel.getSBMLDocument())
    between_growth_test(current_model, config, step="after_draft")
    between_analysis(current_model, config, step="after_draft")

    # gapfilling
    ############
    threshold_add_reacs = config["gapfilling"]["threshold_add_reacs"]
    run_gapfill = False

    # Set up directory for gapfill output, only if any of the gapfillers are activated
    if any([config["gapfilling"]["KEGGapFiller"], config["gapfilling"]["BioCycGapFiller"], config["gapfilling"]["GeneGapFiller"]]):
        current_sub_dir = Path(MISC_DIR, "gapfill")
        try:
            Path(current_sub_dir).mkdir(parents=True, exist_ok=False)
            print(f"Creating new directory {current_sub_dir}")
        except FileExistsError:
            logger.warning(
                f"Given directory {current_sub_dir} already exists. High possibility of files being overwritten."
            )
    
    # KEGGapFiller
    if config["gapfilling"]["KEGGapFiller"]:
        
        logger.info("Filling gaps based on KEGG ...")
        step_start = time.time()
        
        if config["general"]["kegg_organism_id"]:
            run_gapfill = True
            # gapfilling with KEGG via KEGG organism ID
            kgf = KEGGapFiller(config["general"]["kegg_organism_id"])
            kgf.find_missing_genes(current_libmodel)
            kgf.find_missing_reactions(current_model, threshold_add_reacs)
            current_libmodel = kgf.fill_model(
                current_libmodel,
                namespace=config["general"]["namespace"],
                idprefix=config["gapfilling"]["idprefix"],
                formula_check=config["gapfilling"]["formula-check"],
                exclude_dna=config["gapfilling"]["exclude-dna"],
                exclude_rna=config["gapfilling"]["exclude-rna"],
            )
            # save GapFillerReport
            kgf.report(str(current_sub_dir))
            
            # save model
            between_save(
                current_libmodel, MODEL_DIR, "after_KEGG_gapfill", only_modelpath
            )
            current_model = _sbml_to_model(current_libmodel.getSBMLDocument())
        
        elif config["gapfilling"]["KEGGapFiller parameters"]["homologous-id"]:
            logger.info("...using the KEGG ID of a closely related strain.")
            
            run_gapfill = True
            
            # Get BLAST results if they are not passed beforehand
            if config["gapfilling"]["KEGGapFiller parameters"]["blast-result"]:
                blast_result = config["gapfilling"]["KEGGapFiller parameters"]["blast-result"]
            else:
                subprocess.run(
                    [
                        "diamond",
                        "makedb",
                        "--in",
                        config["gapfilling"]["KEGGapFiller parameters"]["homologous-fasta"],
                        "--db",
                        str(Path(MISC_DIR, "gapfill","homologous_strain.dmnd")),
                    ]
                )
                blast_result = run_DIAMOND_blastp(
                    config["general"]["protein_fasta"],
                    db=str(Path(MISC_DIR, "gapfill","homologous_strain.dmnd")),
                    outdir=str(Path(MISC_DIR, "gapfill")),
                    outname="blast_results_homologous_strain.tsv"
                )
            
            # gapfilling with KEGG via alternative KEGG organism ID
            kgf = KEGGapFiller_via_homolog(config["gapfilling"]["KEGGapFiller parameters"]["homologous-id"])
            kgf.find_missing_genes(
                current_libmodel,
                original_fasta=config["general"]["protein_fasta"],
                alternative_fasta=config["gapfilling"]["KEGGapFiller parameters"]["homologous-fasta"],
                blast_result=blast_result,
                cutoff=config["gapfilling"]["KEGGapFiller parameters"]["cutoff"])
            kgf.find_missing_reactions(current_model, threshold_add_reacs)
            current_libmodel = kgf.fill_model(
                current_libmodel,
                namespace=config["general"]["namespace"],
                idprefix=config["gapfilling"]["idprefix"],
                formula_check=config["gapfilling"]["formula-check"],
                exclude_dna=config["gapfilling"]["exclude-dna"],
                exclude_rna=config["gapfilling"]["exclude-rna"],
            )
            # save GapFillerReport
            kgf.report(str(current_sub_dir))
            
            # save model
            between_save(
                current_libmodel, MODEL_DIR, "after_KEGG_gapfill_via_homolog", only_modelpath
            )
            current_model = _sbml_to_model(current_libmodel.getSBMLDocument())
            
        else:
            mes = f"No KEGG organism ID provided. Gapfilling with KEGG will be skipped."
            logger.warning(mes, UserWarning)
        
        step_end = time.time()
        logger.info(f"KEGGapFiller\truntime: {step_end-step_start}s\n")

    # BioCycGapFiller
    if config["gapfilling"]["BioCycGapFiller"]:
        
        logger.info("Filling gaps based on BioCyc ...")
        step_start = time.time()
        
        run_gapfill = True
        
        # Get BLAST results if they are not passed beforehand
        if config["gapfilling"]["BioCycGapFiller parameters"]["homologous-fasta"]:
            homology = True
            
            if config["gapfilling"]["BioCycGapFiller parameters"]["blast-result"]:
                blast_result = config["gapfilling"]["BioCycGapFiller parameters"]["blast-result"]
            else:
                subprocess.run(
                    [
                        "diamond",
                        "makedb",
                        "--in",
                        config["gapfilling"]["BioCycGapFiller parameters"]["homologous-fasta"],
                        "--db",
                        str(Path(MISC_DIR, "gapfill","homologous_strain.dmnd")),
                    ]
                )
                blast_result = run_DIAMOND_blastp(
                    config["general"]["protein_fasta"],
                    db=str(Path(MISC_DIR, "gapfill","homologous_strain.dmnd")),
                    outdir=str(Path(MISC_DIR, "gapfill")),
                    outname="blast_results_homologous_strain.tsv"
                )
        else: 
            blast_result = None
            homology = False
        
        bcgf = BioCycGapFiller(
            config["gapfilling"]["BioCycGapFiller parameters"]["gene-table"],
            config["gapfilling"]["BioCycGapFiller parameters"]["reacs-table"],
            config["gapfilling"]["BioCycGapFiller parameters"]["gff"],
            homology
        )

        bcgf.find_missing_genes(
            current_libmodel,
            fasta_original=config["general"]["protein_fasta"],
            fasta_homolog=config["gapfilling"]["BioCycGapFiller parameters"]["homologous-fasta"],
            blast_result=blast_result,
            cutoff=config["gapfilling"]["BioCycGapFiller parameters"]["cutoff"])
        bcgf.find_missing_reactions(current_model)
        current_libmodel = bcgf.fill_model(
            current_libmodel,
            namespace=config["general"]["namespace"],
            idprefix=config["gapfilling"]["idprefix"],
            formula_check=config["gapfilling"]["formula-check"],
            exclude_dna=config["gapfilling"]["exclude-dna"],
            exclude_rna=config["gapfilling"]["exclude-rna"],
        )
        # save GapFillerReport
        bcgf.report(str(current_sub_dir))
    
        # save model
        if config["gapfilling"]["BioCycGapFiller parameters"]["homologous-fasta"]:
            name = "after_BioCyc_gapfill_via_homolog"
        else: name = "after_BioCyc_gapfill"
        between_save(
            current_libmodel, MODEL_DIR, name, only_modelpath
        )
        current_model = _sbml_to_model(current_libmodel.getSBMLDocument())
        
        step_end = time.time()
        logger.info(f"BioCycGapfiller\truntime: {step_end-step_start}s\n")

    # GeneGapFiller
    if config["gapfilling"]["GeneGapFiller"]:
        
        logger.info("Filling gaps based on genes ...")
        step_start = time.time()
        
        run_gapfill = True
        ggf = GeneGapFiller()
        ggf.find_missing_genes(
            config["gapfilling"]["GeneGapFiller parameters"]["gff"], current_libmodel
        )
        ggf.find_missing_reactions(
            current_model,
            prefix=config["gapfilling"]["idprefix"],
            type_db=config["gapfilling"]["GeneGapFiller parameters"]["type"],
            fasta=config["general"]["protein_fasta"],
            dmnd_db=config["gapfilling"]["GeneGapFiller parameters"]["dmnd-database"],
            map_db=config["gapfilling"]["GeneGapFiller parameters"]["database-mapping"],
            mail=config["tech-resources"]["email"],
            check_NCBI=config["gapfilling"]["GeneGapFiller parameters"]["check-NCBI"],
            threshold_add_reacs=threshold_add_reacs,
            outdir=current_sub_dir,
            sens=config["gapfilling"]["GeneGapFiller parameters"]["sensitivity"],
            cov=config["gapfilling"]["GeneGapFiller parameters"]["coverage"],
            t=config["tech-resources"]["threads"],
            pid=config["gapfilling"]["GeneGapFiller parameters"]["percentage identity"],
        )
        current_libmodel = ggf.fill_model(
            current_libmodel,
            namespace=config["general"]["namespace"],
            idprefix=config["gapfilling"]["idprefix"],
            formula_check=config["gapfilling"]["formula-check"],
            exclude_dna=config["gapfilling"]["exclude-dna"],
            exclude_rna=config["gapfilling"]["exclude-rna"],
        )
        # save GapFillerReport
        ggf.report(str(current_sub_dir))
        
        # save model
        between_save(
            current_libmodel, MODEL_DIR, "after_Gene_gapfill", only_modelpath
        )
        
        step_end = time.time()
        logger.info(f"GeneGapFiller\truntime: {step_end-step_start}s\n")

    current_model = _sbml_to_model(current_libmodel.getSBMLDocument())

    # testing
    if run_gapfill:
        between_growth_test(current_model, config, step="after_gapfill")
        between_analysis(current_model, config, step="after_gapfill")

    # model cleanup - part 1
    ########################

    # duplicates
    # ----------
    
    logger.info("Handling duplicates ...")
    step_start = time.time()
    
    match config["duplicates"]["reactions"]:
        case "remove":
            check_dupl_reac = True
            remove_dupl_reac = True
        case "check":
            check_dupl_reac = True
            remove_dupl_reac = False
        case "skip":
            check_dupl_reac = False
            remove_dupl_reac = False
        case _:
            mes = "Unknown input for duplicates - reactions: will be skipped"
            logger.warning(mes)
            check_dupl_reac = False
            remove_dupl_reac = False

    match config["duplicates"]["metabolites"]:
        case "remove":
            check_dupl_meta = "default"
            remove_dupl_meta = True
        case "check":
            check_dupl_meta = "default"
            remove_dupl_meta = False
        case "skip":
            check_dupl_meta = "skip"
            remove_dupl_meta = False
        case _:
            mes = "Unknown input for duplicates - metabolites: will be skipped"
            logger.warning(mes)
            check_dupl_meta = "skip"
            remove_dupl_meta = False

    current_model = resolve_duplicates(
        current_model,
        check_reac=check_dupl_reac,
        check_meta=check_dupl_meta,
        replace_dupl_meta=remove_dupl_meta,
        remove_unused_meta=config["duplicates"]["remove_unused_metabs"],
        remove_dupl_reac=remove_dupl_reac,
    )
    
    step_end = time.time()
    logger.info(f"Duplicates\truntime: {step_end-step_start}s\n")

    # save model
    between_save(
        current_model, MODEL_DIR, "after_duplicate_removal", only_modelpath
    )

    # in-between testing
    between_growth_test(current_model, config, step="after_duplicate_removal")
    between_analysis(current_model, config, step="after_duplicate_removal")

    current_libmodel = convert_cobra_to_libsbml(current_model, 'notes')

    # ModelPolisher
    ###############
    if config["modelpolisher"]:
        logger.warning('ModelPolisher is currently not maintained and might not work as expected. Use at your own risk.')
        
        logger.info("Running ModelPolisher ...")
        step_start = time.time()
        
        config_mp = {
            "allow-model-to-be-saved-on-server": config["mp"][
                "allow-model-to-be-saved-on-server"
            ],
            "fixing": {"dont-fix": config["mp"]["fixing"]["dont-fix"]},
            "annotation": {
                "bigg": {
                    "annotate-with-bigg": config["mp"]["annotation"]["bigg"][
                        "annotate-with-bigg"
                    ],
                    "include-any-uri": config["mp"]["annotation"]["bigg"][
                        "include-any-uri"
                    ],
                }
            },
        }

        result = run_ModelPolisher(current_libmodel, config_mp)

        if result:
            # Set up directory for result output
            current_sub_dir = Path(MISC_DIR, "modelpolisher", "run_" + result['run_id'])
            try:
                Path(current_sub_dir).mkdir(parents=True, exist_ok=False)
                print(f"Creating new directory {current_sub_dir}")
            except FileExistsError:
                logger.warning(
                    f"Given directory {current_sub_dir} already exists. High possibility of files being overwritten."
                )
            
            if len(result['diff']) > 1:
                pd.DataFrame(result["diff"]).to_csv(
                    Path(current_sub_dir, "diff_mp.csv"),
                    sep=";",
                    header=False,
                )
            else:
                logger.warning(f'{result["diff"]}')
            
            pd.DataFrame(result["pre_validation"]).to_csv(
                Path(current_sub_dir, "pre_validation.csv"),
                sep=";",
                header=True,
            )
            pd.DataFrame(result["post_validation"]).to_csv(
                Path(current_sub_dir, "post_validation.csv"),
                sep=";",
                header=True,
            )

            # save model
            current_libmodel = result["polished_document"].getModel()
            between_save(
                current_libmodel, MODEL_DIR, "after_ModelPolisher", only_modelpath
            )

            current_model = _sbml_to_model(current_libmodel.getSBMLDocument())

            # in-between testing
            between_growth_test(current_model, config, step="after_ModelPolisher")
            between_analysis(current_model, config, step="after_ModelPolisher")
        else:
            logger.warning('No result was produced with ModelPolisher. This step will be skipped.')
            
        step_end = time.time()
        logger.info(f"ModelPolisher\truntime: {step_end-step_start}s\n")

    # Annotations
    #############

    # @BUG/@WARNING: 
    # The correct annotation for pathways gets missing between after the following step and reaction direction.
    # @IDEA Possible solution would be to run polish_annotations & change_all_qualifiers at the end again?
    # KEGGPathwayGroups
    # -----------------
    if config["kegg_pathway_groups"]["add"]:
        
        logger.info("Adding KEGG pathway groups ...")
        step_start = time.time()
        
        missing_list = set_kegg_pathways(
            current_libmodel,
            viaEC=config["kegg_pathway_groups"]["viaEC"],
            viaRC=config["kegg_pathway_groups"]["viaRC"],
        )
        if missing_list:
            with open(
                Path(
                    MISC_DIR,
                    "kegg_pathways",
                    "reac_wo_kegg_pathway_groups.txt",
                ),
                "w",
            ) as outfile:
                for line in missing_list:
                    outfile.write(f"{line}\n")
        # save model
        between_save(
            current_libmodel, MODEL_DIR, "added_KeggPathwayGroups", only_modelpath
        )
        
        step_end = time.time()
        logger.info(f"KEGGPathwayGroups\truntime: {step_end-step_start}s\n")

    # SBOannotator
    # ------------
    logger.info("Running SBOannotator ...")
    step_start = time.time()
    
    current_libmodel = run_SBOannotator(current_libmodel)
    
    step_end = time.time()
    logger.info(f"SBOannotator\truntime: {step_end-step_start}s\n")
    
    between_save(
        current_libmodel, MODEL_DIR, "SBOannotated", only_modelpath
    )
    
    # ECO Terms
    # ---------
    
    logger.info("Adding ECO terms ...")
    step_start = time.time()
    
    current_libmodel = add_ECO_terms(current_libmodel)
    change_all_qualifiers(current_libmodel, lab_strain=config["cm-polish"]["is_lab_strain"])

    step_end = time.time()
    logger.info(f"\truntime: {step_end-step_start}s\n")

    between_save(current_libmodel, MODEL_DIR, "added_ECO_terms", only_modelpath)

    current_model = _sbml_to_model(current_libmodel.getSBMLDocument())
    between_analysis(current_model, config, step="after_annotation")
    
    # Model cleanup - part 2
    ########################

    # reaction direction
    # ------------------
    if config["reaction_direction"]:
        logger.info("Checking reaction direction with BioCyc ...")
        step_start = time.time()
        
        current_model = check_direction(current_model, config["reaction_direction"])

        step_end = time.time()
        logger.info(f"Reaction direction\truntime: {step_end-step_start}s\n")

    # save model
    between_save(
        current_model, MODEL_DIR, "after_reac_direction_change", only_modelpath
    )
    
    # MCC
    # ---
    logger.info("Running MCC ...")
    step_start = time.time()
    
    current_model = perform_mcc(
        current_model, Path(MISC_DIR, "mcc"), apply=True
    )
    
    step_end = time.time()
    logger.info(f"MCC\truntime: {step_end-step_start}s\n")

    # find and solve energy generating cycles
    # ---------------------------------------
    
    logger.info("Finding (and resolving) EGCs ...")
    step_start = time.time()
    
    results = None
    match config["EGCs"]["solver"]:
        # greedy solver
        case "greedy":
            logger.info("Using GreedyEGCSolver...")
            solver = egcs.GreedyEGCSolver()
            # automatically uses c,e as compartments
            results = solver.solve_egcs(
                current_model, namespace=config["general"]["namespace"]
            )
            if results:
                logger.info("results:")
                for k, v in results.items():
                    logger.info(f"\t{k}: {v}")

        # no solver = EGCs will only be reported
        case _:
            solver = egcs.EGCSolver()
            logger.info(f"\tFound EGCs:\n")
            # automatically uses c,e as compartments
            logger.info(
                f'\t{solver.find_egcs(current_model,with_reacs=True,namespace=config["general"]["namespace"])}'
            )
            
    step_end = time.time()
    logger.info(f"EGCs\truntime: {step_end-step_start}s\n")

    if results:
        between_save(
            current_model, MODEL_DIR, "after_egc_fix", only_modelpath
        )
        current_libmodel = convert_cobra_to_libsbml(current_model, 'notes')
        # in-between testing
        between_growth_test(current_model, config, step="after_egcs")
        between_analysis(current_model, config, step="after_egcs")

    # BOF
    # ---
    
    logger.info("Improve biomass objective function ...")
    step_start = time.time()
    
    # BOFdat - optional
    if config["BOF"]["run_bofdat"]:
        check_bof = test_biomass_presence(current_model)
        if check_bof:
            current_model.reactions.get_by_id(check_bof[0]).reaction = adjust_BOF(
                config["BOF"]["bofdat_params"]["full_genome_sequence"],
                current_model,
                dna_weight_fraction=config["BOF"]["bofdat_params"][
                    "dna_weight_fraction"
                ],
                weight_frac=config["BOF"]["bofdat_params"]["weight_fraction"],
            )
        else:
            bof_reac = Reaction("Biomass_BOFdat")
            bof_reac.name = "Biomass objective function created by BOFdat"
            current_model.add_reactions([bof_reac])
            current_model.reactions.get_by_id(bof_reac).reaction = adjust_BOF(
                config["BOF"]["bofdat_params"]["full_genome_sequence"],
                current_model,
                dna_weight_fraction=config["BOF"]["bofdat_params"][
                    "dna_weight_fraction"
                ],
                weight_frac=config["BOF"]["bofdat_params"]["weight_fraction"],
            )
    # check and normalise
    current_model = check_normalise_biomass(current_model)
    
    step_end = time.time()
    logger.info(f"BOF\truntime: {step_end-step_start}s\n")

    # save model
    between_save(current_model, MODEL_DIR, "after_BOF", only_modelpath)
    
    # prune model 
    # -----------
    # (remove mass unbalanced reactions to ensure model consistency)
    prune_mass_unbalanced_reacs(current_model)

    # Clean-up & save the final model
    # Set up directory for final result
    current_sub_dir = Path(MISC_DIR, "clean-up")
    try:
        Path(current_sub_dir).mkdir(parents=True, exist_ok=False)
        print(f"Creating new directory {current_sub_dir}")
    except FileExistsError:
        logger.warning(
            f"Given directory {current_sub_dir} already exists. High possibility of files being overwritten."
        )
    current_libmodel = convert_cobra_to_libsbml(current_model, 'notes')
    # move qualifiers to correct biological or model qualifier type 
    current_libmodel = polish_annotations(current_libmodel, True, str(current_sub_dir))
    current_libmodel = change_all_qualifiers(current_libmodel, lab_strain=config["cm-polish"]["is_lab_strain"])
    between_save(current_libmodel, MODEL_DIR, "final_model", only_modelpath)

    # analysis
    ##########
    
    logger.info("Analysing the final model ...")
    step_start = time.time()

    # stats
    # -----
    logger.info("\tGenerating statistics report ...")
    stats_report = ModelInfoReport(current_model)
    current_sub_dir = Path(MISC_DIR, "stats", "final")
    try:
        Path(current_sub_dir).mkdir(parents=True, exist_ok=False)
    except FileExistsError:
        logger.warning(
            f"Given directory {current_sub_dir} already exists. High possibility of files being overwritten."
        )
    stats_report.save(
        Path(MISC_DIR, "stats", "final"),
        color_palette=config["general"]["colours"],
    )

    # kegg pathway
    # ------------
    logger.info("\tGenerating KEGG pathway report ...")
    pathway_report = kegg_pathway_analysis(current_model)
    pathway_report.save(
        Path(MISC_DIR, "kegg_pathways"),
        colors=config["general"]["colours"],
    )

    # sbo term
    # --------
    logger.info("\tAnalysing SBO terms ...")
    sboreport = SBOTermReport(current_libmodel)
    sboreport.save(str(Path(MISC_DIR)))

    # memote
    # ------
    logger.info("\tGenerating memote report ...")
    run_memote(
        current_model,
        "html",
        save_res=Path(MISC_DIR, "memote", "final_memote.html"),
    )

    # growth
    # ------
    logger.info("\t Performing final growth tests...")
    between_growth_test(current_model, config, step="final")

    # auxotrophies
    # ------------
    logger.info("\tChecking for amino acids auxotrophies ...")
    media_list, suppl_list = load_media(config["input"]["mediapath"])
    auxo_report = growth.test_auxotrophies(
        current_model, media_list, suppl_list, config["general"]["namespace"]
    )
    auxo_report.save(
        Path(MISC_DIR, "auxotrophy"),
        color_palette=config["general"]["colours"],
    )
    
    step_end = time.time()
    logger.info(f"Final analysis\truntime: {step_end-step_start}s\n")

    total_time_e = time.time()
    logger.info(f'Total run time: {total_time_e - total_time_s}')


# run for multiple models
@implement
def wrapper():
    """Run given settings for the CarveMe-ModelPolisher-based (CMPB) workflow on multiple models.
    
    .. note::
    
        This function is not yet implemented. It is a placeholder for future development.
    """
    # dev notes:
    # - Add option to have specific colour list per model for plots
    # for the comparison when runnin this on multiple models
    # - Maybe get models at first and then add model IDs to every save filename?
    pass
