#!/usr/bin/env python
"""Workflow call for the CarveMe-ModelPolisher-based (CMPB) workflow.
"""

__author__ = "Tobias Fehrenbach, Famke Baeuerle, Gwendolyn O. Döbel and Carolin Brune"

################################################################################
# requirements
################################################################################

from datetime import date
import logging
import model_polisher as mp
import os
import pandas as pd
from pathlib import Path
from typing import Union

import warnings

from cobra import Reaction,Model
from libsbml import Model as libModel

from refinegems.analysis import growth
from refinegems.analysis.investigate import plot_rea_sbo_single
from refinegems.classes.reports import ModelInfoReport
from refinegems.utility.util import test_biomass_presence
from refinegems.curation.biomass import check_normalise_biomass
from refinegems.curation.charges import correct_charges_modelseed
from refinegems.curation.curate import resolve_duplicates
from refinegems.classes.gapfill import KEGGapFiller, BioCycGapFiller, GeneGapFiller
from refinegems.curation.pathways import kegg_pathways, kegg_pathway_analysis
from refinegems.curation.polish import polish, check_direction
from refinegems.utility.entities import resolve_compartment_names, are_compartment_names_valid
from refinegems.utility.connections import run_memote, perform_mcc, adjust_BOF, run_SBOannotator
from refinegems.utility.io import load_model, write_model_to_file
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

    # ....................................................
    # @TODO / @IDEA
    # use temp folder or report all model/in-between steps
    # what to write in the log file
    # e.g. runtimes, warnings, hints and more
    # ....................................................

# dev notes
#   in the run function: current_model means the cobrapy model, 
#   while current_libmodel means the libsbml model
# @TEST
def run(configpath:Union[str,None]=None):
    """Run the CarveMe-ModelPolisher-based (CMPB) workflow.

    Args:
        - configpath (Union[str,None]): 
            The Path to a configuration file. If none given, prompts the user to enter the needed 
            information step by step.
            Defaults to None.
    """

    def between_save(model:Union[Model,libModel], dir:Union[None,str]=None, 
                     label:Union[None,str]=None, 
                     only_modelpath:Union[None,str]=None) -> str:
        """Helper function to save the model in between the different steps of the workflow

        Args:
            - model (Union[Model,libModel]): 
                The current model, either loaded with libsbml or COBRApy.
            - dir (Union[None,str], optional): 
                Output directory. Defaults to None.
            - label (Union[None,str], optional): 
                Name of the step. Defaults to None.
            - only_modelpath (Union[None,str], optional): 
                Path to save the model to, if only one model will be saved. 
                If set, ignores dir and label.
                Defaults to None.

        Returns:
            str: 
                The current model path.
        """

        name = model.id if isinstance(model, Model) else model.getId()
        if not only_modelpath:
            write_model_to_file(model, str(Path(dir,'cmpb_out','models',f'{name}_{label}.xml')))     
            return Path(dir,'cmpb_out','models',f'{name}_{label}.xml')
        else:
            write_model_to_file(model, str(only_modelpath))
            return only_modelpath
        

    def between_growth_test(model: Model, cfg:dict, step:str):
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
            growth_report = growth.growth_analysis(model, cfg['input']['mediapath'], 
                                                namespace=cfg['general']['namespace'], 
                                                retrieve='report')
            growth_report.save(Path(cfg['general']["dir"], 'cmpb_out', 'misc', 'growth', step),
                               color_palette=config['general']['colours']) 
        else:
            mes = f'No growth/biomass function detected, growth simulation for step {step} will be skipped.'
            warnings.warn(mes)


    def between_analysis(model: Model, cfg:dict, step:str):
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
        if cfg['general']['memote_always_on']:
            run_memote(model, 'html', 
                    return_res=False, 
                    save_res=Path(cfg['general']["dir"], 'cmpb_out', 'misc','memote',f'{step}.html'),
                    verbose=False)
        if cfg['general']['stats_always_on']:
            report = ModelInfoReport(model)
            Path(dir,"cmpb_out",'misc', 'stats', step).mkdir(parents=True, exist_ok=False)
            report.save(Path(cfg['general']["dir"],'cmpb_out', 'misc', 'stats', step),
                        color_palette=config['general']['colours']) 


    # setup phase
    #############

    # load config
    # -----------
    if not configpath:
        config = save_cmpb_user_input()
    else:
        config = validate_config(configpath, 'cmpb')

    # set up model name
    # -----------------
    if config['general']['authorinitials'] is not None and config['general']['organism'] is not None and config['general']['strainid'] is not None:
        modelname = 'i'+config['general']['organism']+str(config['general']['strainid'])+config['general']['authorinitials']+str(date.today().year).removeprefix('20')
    elif config['general']['modelname'] is not None:
        modelname = config['general']['modelname']
    else:
        print("No values given for the standard name for a model. Default name will be used.")
        modelname = "model_"+str(date.today().year).removeprefix('20')

    dir = config['general']['dir']
    only_modelpath = None
    if not config['general']['save_all_models']:
        only_modelpath = Path(dir,'cmpb_out',f'{modelname}.xml') 

    # create directory structure
    # --------------------------
    build_data_directories('cmpb',dir)
    
    # create log
    # ----------
    today = date.today().strftime("%Y%m%d")
    log_file = Path(dir, 'cmpb_out', 'logs', f'specimen_cmpb_{str(today)}.log')
    handler = logging.handlers.RotatingFileHandler(log_file, 
                                                    mode='w', 
                                                    backupCount=10, 
                                                    encoding='utf-8', 
                                                    delay=0)
    handler.setFormatter(logging.Formatter("{levelname} \t {name} \t {message}", 
                                           style="{",))
    logger.addHandler(handler)

    # CarveMe
    #########
    if not config['input']['modelpath']:
        if config['carveme']['gram'] == "grampos" or config['carveme']['gram'] == "gramneg":
            os.system(f"carve {config['general']['protein_fasta']} --solver scip -u {config['carveme']['gram']} -o {dir}\cmpb_out\models\{modelname}.xml")
        else: 
            os.system(f"carve {config['general']['protein_fasta']} --solver scip -o {dir}\cmpb_out\models\{modelname}.xml")
        config['input']['modelpath'] = dir+f'\cmpb_out\models\{modelname}.xml'
    current_modelpath = config['input']['modelpath']
    
    
    # validate compartments of the loaded model
    ###########################################
    current_model = load_model(str(current_modelpath),'cobra')
    # check, if compartment names are valid
    if not are_compartment_names_valid(current_model):
        # if not, attempt to fix them
        resolve_compartment_names(current_model)
        # save validated model
        current_modelpath = between_save(current_model, dir, 'validated_comps',only_modelpath)
            
    # reload into libsbml
    current_libmodel = load_model(str(current_modelpath),'libsbml')    
    

    # CarveMe correction
    ####################
    
    # check, if input is a CarveMe model
    if 'CarveMe' in current_libmodel.getAnnotationString():
        Path(dir,"cmpb_out",'misc', 'wrong_annotations').mkdir(parents=True, exist_ok=False)
        current_libmodel = polish(current_libmodel, 
                                  email = config['tech-resources']['email'],
                                  id_db = config['general']['namespace'],
                                  gff = config['general']['gff'],
                                  protein_fasta = config['general']['protein_fasta'],
                                  lab_strain = config['cm-polish']['is_lab_strain'],
                                  kegg_organism_id = config['general']['kegg_organism_id'],
                                  path = Path(dir,'cmpb_out','misc','wrong_annotations'))
        # rg correct charges
        current_libmodel, mult_charges_dict = correct_charges_modelseed(current_libmodel)
        mult_charges_tab = pd.DataFrame.from_dict(mult_charges_dict, orient='index')
        mult_charges_tab.to_csv(Path(dir,'cmpb_out','misc','reac_with_mult_charges.tsv'), sep='\t')
        
        # save model
        current_modelpath = between_save(current_libmodel, dir, 'after_CarveMe_correction',only_modelpath)
        

    # growth test
    # -----------
    current_model = load_model(str(current_modelpath),'cobra')
    between_growth_test(current_model,config,step='after_draft')
    between_analysis(current_model, config, step='after_draft')

    # gapfilling
    ############
    threshold_add_reacs = config['gapfilling']['threshold_add_reacs']
    run_gapfill = False
    # KEGGapFiller
    if config['gapfilling']['KEGGapFiller']:
        if config['general']['kegg_organism_id']:
            run_gapfill = True
            # gapfilling with KEGG via KEGG organism ID
            kgf = KEGGapFiller(config['general']['kegg_organism_id'])
            kgf.find_missing_genes(current_libmodel)
            kgf.find_missing_reactions(current_model, threshold_add_reacs)
            current_libmodel = kgf.fill_model(current_libmodel, 
                                              namespace = config['general']['namespace'],
                                              idprefix = config['gapfilling']['idprefix'],
                                              formula_check = config['gapfilling']['formula-check'],
                                              exclude_dna = config['gapfilling']['formula-check'],
                                              exclude_rna = config['gapfilling']['exclude-rna']
                                              )
            # save model
            current_modelpath = between_save(current_libmodel, dir, 'after_KEGG_gapfill',only_modelpath)
            current_model = load_model(str(current_modelpath),'cobra')
            
        else:
            mes = f'No KEGG organism ID provided. Gapfilling with KEGG will be skipped.'
            raise warnings.warn(mes,UserWarning)
    
    # BioCycGapFiller
    if config['gapfilling']['BioCycGapFiller']:
        run_gapfill = True
        bcgf = BioCycGapFiller(config['gapfilling']['BioCycGapFiller parameters']['gene-table'],
                               config['gapfilling']['BioCycGapFiller parameters']['reacs-table'],
                               config['gapfilling']['BioCycGapFiller parameters']['gff']
                               )

        bcgf.find_missing_genes(current_libmodel)
        bcgf.find_missing_reactions(current_model, threshold_add_reacs)
        current_libmodel = kgf.fill_model(current_libmodel, 
                                          namespace = config['general']['namespace'],
                                          idprefix = config['gapfilling']['idprefix'],
                                          formula_check = config['gapfilling']['formula-check'],
                                          exclude_dna = config['gapfilling']['formula-check'],
                                          exclude_rna = config['gapfilling']['exclude-rna']
                                         )
        # save model
        current_modelpath = between_save(current_libmodel, dir, 'after_BioCyc_gapfill',only_modelpath)
        current_model = load_model(str(current_modelpath),'cobra')
        
    # GeneGapFiller
    if config['gapfilling']['GeneGapFiller']:
        run_gapfill = True
        ggf = GeneGapFiller()
        ggf.find_missing_genes(config['gapfilling']['GeneGapFiller parameters']['gff'],
                               current_libmodel)
        ggf.find_missing_reactions(current_model,
                                   prefix=config['gapfilling']['idprefix'],
                                   type_db=config['gapfilling']['GeneGapFiller parameters']['type'],
                                   fasta = config['general']['protein_fasta'],
                                   dmnd_db = config['gapfilling']['GeneGapFiller parameters']['dmnd-database'],
                                   map_db = config['gapfilling']['GeneGapFiller parameters']['database-mapping'],
                                   mail = config['tech-resources']['email'],
                                   check_NCBI = config['gapfilling']['GeneGapFiller parameters']['check-NCBI'],
                                   threshold_add_reacs = threshold_add_reacs,
                                   outdir = Path(dir,"cmpb_out",'misc', 'gapfill'), # @DISCUSSION or would a subfolder be better?
                                   sens = config['gapfilling']['GeneGapFiller parameters']['sensitivity'],
                                   cov = config['gapfilling']['GeneGapFiller parameters']['coverage'],
                                   t = config['tech-resources']['threads'],
                                   pid = config['gapfilling']['GeneGapFiller parameters']['percentage identity'] 
                               )
        current_libmodel = ggf.fill_model(current_libmodel, 
                                          namespace = config['general']['namespace'],
                                          idprefix = config['gapfilling']['idprefix'],
                                          formula_check = config['gapfilling']['formula-check'],
                                          exclude_dna = config['gapfilling']['formula-check'],
                                          exclude_rna = config['gapfilling']['exclude-rna']
                                         )
        # save model            
        current_modelpath = between_save(current_libmodel, dir, 'after_Gene_gapfill',only_modelpath)
        
    current_model = load_model(str(current_modelpath),'cobra')

    # testing
    if run_gapfill:
        between_growth_test(current_model,config,step='after_gapfill')
        between_analysis(current_model, config, step='after_gapfill')
        

    # ModelPolisher
    ###############
    if config['modelpolisher']:
        config_mp = {"allow-model-to-be-saved-on-server": config["mp"]["allow-model-to-be-saved-on-server"], 
                         "fixing": {"dont-fix": config["mp"]["fixing"]["dont-fix"]},
                         "annotation": {"bigg": {"annotate-with-bigg": config["mp"]["annotation"]["bigg"]["annotate-with-bigg"], 
                                                 "include-any-uri": config["mp"]["annotation"]["bigg"]["include-any-uri"]}}}
        
        result = mp.polish_model_file(current_modelpath, config_mp)
        
        # @DISCUSSION Should the run-id be saved somewhere for debugging purposes? result['run_id']
        pd.DataFrame(result['diff']).to_csv(Path(dir,'cmpb_out','misc','modelpolisher','diff_mp.csv'), sep=';', header=False)
        pd.DataFrame(result['pre_validation']).to_csv(Path(dir,'cmpb_out','misc','modelpolisher','pre_validation.csv'), sep=';', header=True)
        pd.DataFrame(result['post_validation']).to_csv(Path(dir,'cmpb_out','misc','modelpolisher','post_validation.csv'), sep=';', header=True)

        # save model
        current_modelpath = between_save(current_libmodel, dir, 'after_ModelPolisher',only_modelpath)

        current_model = load_model(str(current_modelpath),'cobra')

        # in-between testing
        between_growth_test(current_model,config,step='after_ModelPolisher')
        between_analysis(current_model, config, step='after_ModelPolisher')

    # Annotations
    #############

    # KEGGPathwayGroups
    # -----------------
    if config['kegg_pathway_groups']:
        current_libmodel, missing_list = kegg_pathways(current_modelpath)
        with open(Path(dir, 'cmpb_out', 'misc', 'kegg_pathway', 'reac_wo_kegg_pathway_groups.txt'), 'w') as outfile:
            for line in missing_list:
                outfile.write(f'{line}\n')
        # save model
        current_modelpath = between_save(current_libmodel, dir, 'added_KeggPathwayGroups',only_modelpath)


    # SBOannotator
    # ------------
    current_libmodel = load_model(str(current_modelpath),'libsbml')
    current_libmodel = run_SBOannotator(current_libmodel)
    current_modelpath = between_save(current_libmodel, dir, 'SBOannotated',only_modelpath)


    current_model = load_model(str(current_modelpath),'cobra')
    between_analysis(current_model,config,step='after_annotation')
    
    # model cleanup
    ###############
    
    # duplicates
    # ----------
    match config['duplicates']['reactions']:
        case 'remove':
            check_dupl_reac = True
            remove_dupl_reac = True
        case 'check':
            check_dupl_reac = True
            remove_dupl_reac = False
        case 'skip':
            check_dupl_reac = False
            remove_dupl_reac = False
        case _:
            mes = 'Unknown input for duplicates - reactions: will be skipped'
            warnings.warn(mes)
            check_dupl_reac = False
            remove_dupl_reac = False
        
    match config['duplicates']['metabolites']:
        case 'remove':
            check_dupl_meta = 'default'
            remove_dupl_meta = True
        case 'check':
            check_dupl_meta = 'default'
            remove_dupl_meta = False
        case 'skip':
            check_dupl_meta = 'skip'
            remove_dupl_meta = False
        case _:
            mes = 'Unknown input for duplicates - metabolites: will be skipped'
            warnings.warn(mes)
            check_dupl_meta = 'skip'
            remove_dupl_meta = False

    current_model = resolve_duplicates(current_model, check_reac=check_dupl_reac, 
                       check_meta=check_dupl_meta, 
                       replace_dupl_meta=remove_dupl_meta, 
                       remove_unused_meta=config['duplicates']['remove_unused_metabs'], 
                       remove_dupl_reac=remove_dupl_reac)

    # save model
    current_modelpath = between_save(current_model, dir, 'after_duplicate_removal',only_modelpath)
    
    # in-between testing
    between_growth_test(current_model,config,step='after_duplicate_removal')
    between_analysis(current_model,config,step='after_duplicate_removal')

    # reaction direction
    # ------------------
    # @TEST enabling this 
    if config['reaction_direction']:
        current_model = load_model(str(current_modelpath),'cobra')
        current_model = check_direction(current_model, config['reaction_direction'])

    # save model
    current_modelpath = between_save(current_model, dir, 'after_reac_direction_change',only_modelpath)
    
    # find and solve energy generating cycles
    # ---------------------------------------
    current_model = load_model(str(current_modelpath),'cobra')
    results = None
    match config['EGCs']['solver']:
        # greedy solver
        case 'greedy':
            print('Using GreedyEGCSolver...')
            solver = egcs.GreedyEGCSolver()
            results = solver.solve_egcs(current_model,namespace=config['general']['namespace']) # @NOTE automatically uses c,p as compartments 
            if results: 
                logger.info('results:')
                for k,v in results.items():
                    logger.info(f'\t{k}: {v}')
        
        # no solver = EGCs will only be reported
        case _:
            solver = egcs.EGCSolver()
            logger.info(f'\tFound EGCs:\n')
            logger.info(f'\t{solver.find_egcs(current_model,with_reacs=True,namespace=config["general"]["namespace"])}') # @NOTE automatically uses c,p as compartments 
       
    if results:
        current_modelpath = between_save(current_model, dir, 'after_egc_fix',only_modelpath)
        current_libmodel = load_model(str(current_modelpath),'libsbml')
        # in-between testing
        between_growth_test(current_model,config,step='after_egcs')
        between_analysis(current_model,config,step='after_egcs')
            
    # BOF
    # ---
    # BOFdat - optional
    if config['BOF']['run_bofdat']:
        check_bof = test_biomass_presence(current_model)
        if check_bof:
            current_model.reactions.get_by_id(check_bof[0]).reaction = adjust_BOF(config['BOF']['bofdat_params']['full_genome_sequence'], 
                                str(current_modelpath),
                                current_model, 
                                dna_weight_fraction = config['BOF']['bofdat_params']['dna_weight_fraction'], 
                                weight_frac = config['BOF']['bofdat_params']['weight_fraction']) 
            current_model = check_normalise_biomass(current_model)
        else: 
            bof_reac = Reaction('Biomass_BOFdat')
            bof_reac.name = 'Biomass objective function created by BOFdat'
            current_model.add_reactions([bof_reac])
            current_model.reactions.get_by_id(bof_reac).reaction = adjust_BOF(config['BOF']['bofdat_params']['full_genome_sequence'], 
                                str(current_modelpath),
                                current_model, 
                                dna_weight_fraction = config['BOF']['bofdat_params']['dna_weight_fraction'], 
                                weight_frac = config['BOF']['bofdat_params']['weight_fraction'])
            current_model = check_normalise_biomass(current_model)
    # check and normalise 
    else:
        current_model = check_normalise_biomass(current_model)

    # save model
    current_modelpath = between_save(current_model, dir, 'after_BOF',only_modelpath)
    
    # MCC
    # ---
    current_model = perform_mcc(current_model, Path(dir,'cmpb_out','misc','mcc'),apply=True)

    # save the final model
    current_modelpath = between_save(current_model, dir, 'final_model',only_modelpath)

    # analysis
    ##########
    current_model = load_model(str(current_modelpath),'cobra')
    current_libmodel = load_model(str(current_modelpath),'libsbml')

    # stats
    # -----
    stats_report = ModelInfoReport(current_model)
    Path(dir,'cmpb_out', 'misc','stats','final').mkdir(parents=True, exist_ok=False)
    stats_report.save(Path(dir,'cmpb_out', 'misc','stats','final'),
                      color_palette=config['general']['colours']) 
    
    # kegg pathway
    # ------------
    pathway_report = kegg_pathway_analysis(current_model)
    pathway_report.save(Path(dir,'cmpb_out','misc','kegg_pathway'), 
                        colors=config['general']['colours']) 
    
    # sbo term
    # --------
    fig = plot_rea_sbo_single(current_libmodel)
    fig.savefig(Path(dir,'cmpb_out','misc','sbo_gene.png'), dpi=400)

    # memote
    # ------
    run_memote(current_model, 'html', save_res=Path(dir, 'cmpb_out','misc','memote','final_memote.html'))

    # growth
    # ------
    between_growth_test(current_model,config,step='final')

    # auxotrophies
    # ------------
    media_list = growth.read_media_config(config['input']['mediapath'])
    auxo_report = growth.test_auxotrophies(current_model, media_list[0], media_list[1], config['general']['namespace'])
    auxo_report.save(Path(dir,'cmpb_out','misc','auxotrophy'),
                     color_palette=config['general']['colours']) 


####### IDEAS below ####

# @TODO 
# - Add option to have specific colour list per model for plots
    # for the comparison when runnin this on multiple models
# - Maybe get models at first and then add model IDs to every save filename?
# - Add optional FROG report at end of pipeline

# run for multiple models
@implement
def wrapper():
    """Run given settings for the CarveMe-ModelPolisher-based (CMPB) workflow on multiple models.
    """
    pass