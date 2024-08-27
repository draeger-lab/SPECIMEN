#!/usr/bin/env python
"""Workflow call for the CarveMe-ModelPolisher-based (CMPB) workflow.
"""

__author__ = "Tobias Fehrenbach, Famke Baeuerle, Gwendolyn O. Döbel and Carolin Brune"

################################################################################
# requirements
################################################################################

import logging
import pandas as pd
from datetime import date
from pathlib import Path
from typing import Union

import warnings
import yaml

from cobra import Reaction,Model
from libsbml import readSBML
import subprocess

from refinegems.analysis import growth
from refinegems.analysis.investigate import plot_rea_sbo_single
from refinegems.classes.reports import ModelInfoReport
from refinegems.utility.util import test_biomass_presence
from refinegems.curation.biomass import check_normalise_biomass
from refinegems.curation.charges import correct_charges_modelseed
from refinegems.curation.curate import resolve_duplicates
from refinegems.classes.gapfill import KEGGapFiller, BioCycGapFiller, GeneGapFiller
from refinegems.curation.pathways import kegg_pathways, kegg_pathway_analysis
from refinegems.curation.polish import polish
from refinegems.utility.connections import run_memote, perform_mcc, adjust_BOF, run_SBOannotator
from refinegems.utility.io import load_model, write_model_to_file
from refinegems.developement.decorators import implement

from ..util.set_up import save_cmpb_user_input, validate_config

################################################################################
# functions
################################################################################

    # ....................................................
    # @TODO / @IDEAS
    # use temp folder or report all model/in-between steps
    # what to write in the log file
    # ....................................................

# dev notes
#   in the run function: current_model means the cobrapy model, 
#   while current_libmodel means the libsbml model

def run(configpath:Union[str,None]=None):
    """Run the CarveMe-ModelPolisher-based (CMPB) workflow.

    Args:
        - configpath (Union[str,None]): 
            The Path to a configuration file. If none given, prompts the user to enter the needed 
            information step by step.
            Defaults to None.
    """

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

    dir = config['general']['dir']
    if not config['general']['save_all_models']:
        only_modelpath = Path(dir,'cmpb_out','model.xml') 

    # create directory structure
    # --------------------------

    Path(dir,"cmpb_out").mkdir(parents=True, exist_ok=False)                          # cmpb_out
    Path(dir,"cmpb_out",'models').mkdir(parents=True, exist_ok=False)                 #   |- models
    Path(dir,"cmpb_out",'logs').mkdir(parents=True, exist_ok=False)                   #   |- logs
    Path(dir,"cmpb_out",'misc').mkdir(parents=True, exist_ok=False)                   #   |- misc
    Path(dir,"cmpb_out",'misc', 'memote').mkdir(parents=True, exist_ok=False)         #      |- memote
    Path(dir,"cmpb_out",'misc', 'mcc').mkdir(parents=True, exist_ok=False)            #      |- mcc
    Path(dir,"cmpb_out",'misc', 'gapfill').mkdir(parents=True, exist_ok=False)        #      |- gapfill
    Path(dir,"cmpb_out",'misc', 'growth').mkdir(parents=True, exist_ok=False)         #      |- growth
    Path(dir,"cmpb_out",'misc', 'stats').mkdir(parents=True, exist_ok=False)          #      |- stats
    Path(dir,"cmpb_out",'misc', 'kegg_pathway').mkdir(parents=True, exist_ok=False)   #      |- kegg_pathways
    Path(dir,"cmpb_out",'misc', 'auxotrophy').mkdir(parents=True, exist_ok=False)     #      |- auxothrophy

    # create log
    # ----------
    today = date.today().strftime("%Y%m%d")
    log_file = Path(dir, 'cmpb_out', 'logs', f'specimen_cmpb_{str(today)}.log')

    # CarveMe
    #########
    if not config['input']['modelpath']:
        if config['carveme']['gram'] == "grampos" or config['carveme']['gram'] == "gramneg":
            subprocess.run(["carve", config['general']['protein_fasta'], "--solver", "scip", '-u', config['carveme']['gram'], "-o", dir+"\cmpb_out\models\Draft.xml"])
        else: 
            subprocess.run(["carve", config['general']['protein_fasta'], "--solver", "scip", "-o", dir+"\cmpb_out\models\Draft.xml"])
        config['input']['modelpath'] = dir+'\cmpb_out\models\Draft.xml'
    current_modelpath = config['input']['modelpath']

    # CarveMe correction
    ####################
    current_libmodel = load_model(str(current_modelpath),'libsbml')
    # check, if input is a CarveMe model
    if 'CarveMe' in current_libmodel.getAnnotationString():
        Path(dir,"cmpb_out",'misc', 'wrong_annotations').mkdir(parents=True, exist_ok=False)
        current_libmodel = polish(current_libmodel, 
                                  email = config['tech-resources']['email'],
                                  id_db = config['general']['namespace'],
                                  refseq_gff = config['general']['refseq_gff'],
                                  protein_fasta = config['general']['protein_fasta'],
                                  lab_strain = config['cm-polish']['is_lab_strain'],
                                  kegg_organism_id = config['general']['kegg_organism_id'],
                                  path = Path(dir,'cmpb_out','misc','wrong_annotations'))
        # rg correct charges
        current_libmodel, mult_charges_dict = correct_charges_modelseed(current_libmodel)
        mult_charges_tab = pd.DataFrame.from_dict(mult_charges_dict, orient='index')
        mult_charges_tab.to_csv(Path(dir,'cmpb_out','misc','reac_with_mult_charges.tsv'), sep='\t')
        
        # save model
        if config['general']['save_all_models']:
            write_model_to_file(current_libmodel, str(Path(dir,'cmpb_out','models',f'{current_libmodel.getId()}_after_CarveMe_correction.xml')))     
            current_modelpath = Path(dir,'cmpb_out','models',f'{current_libmodel.getId()}_after_CarveMe_correction.xml')
        else:
            write_model_to_file(current_libmodel, str(only_modelpath))
            current_modelpath = only_modelpath

    # growth test
    # -----------
    current_model = load_model(str(current_modelpath),'cobra')
    between_growth_test(current_model,config,step='after_draft')
    between_analysis(current_model, config, step='after_draft')


    # gapfilling
    ############
    run_gapfill = False
    # KEGGapFiller
    if config['gapfilling']['KEGGapFiller']:
        if config['general']['kegg_organism_id']:
            run_gapfill = True
            # gapfilling with KEGG via KEGG organism ID
            kgf = KEGGapFiller(config['general']['kegg_organism_id'])
            kgf.find_missing_genes(current_libmodel)
            kgf.find_missing_reactions(current_model)
            current_libmodel = kgf.fill_model(current_libmodel, 
                                              namespace = config['general']['namespace'],
                                              idprefix = config['gapfilling']['idprefix'],
                                              formula_check = config['gapfilling']['formula-check'],
                                              exclude_dna = config['gapfilling']['formula-check'],
                                              exclude_rna = config['gapfilling']['exclude-rna']
                                              )
            # save model
            if config['general']['save_all_models']:
                write_model_to_file(current_libmodel, str(Path(dir,'cmpb_out','models',f'{current_libmodel.getId()}_after_KEGG_gapfill.xml')))     
                current_modelpath = Path(dir,'cmpb_out','models',f'{current_libmodel.getId()}_after_KEGG_gapfill.xml')
            else:
                write_model_to_file(current_libmodel, str(only_modelpath))
            current_model = load_model(current_modelpath,'cobra')
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
        # @TEST if it works
        bcgf.find_missing_genes(current_libmodel)
        bcgf.find_missing_reactions(current_model)
        current_libmodel = kgf.fill_model(current_libmodel, 
                                          namespace = config['general']['namespace'],
                                          idprefix = config['gapfilling']['idprefix'],
                                          formula_check = config['gapfilling']['formula-check'],
                                          exclude_dna = config['gapfilling']['formula-check'],
                                          exclude_rna = config['gapfilling']['exclude-rna']
                                         )
        # save model
        if config['general']['save_all_models']:
            write_model_to_file(current_libmodel, str(Path(dir,'cmpb_out','models',f'{current_libmodel.getId()}_after_BioCyc_gapfill.xml')))     
            current_modelpath = Path(dir,'cmpb_out','models',f'{current_libmodel.getId()}_after_BioCyc_gapfill.xml')
        else:
            write_model_to_file(current_libmodel, str(only_modelpath))
        current_model = load_model(current_modelpath,'cobra')
        
    # GeneGapFiller
    if config['gapfilling']['GeneGapFiller']:
        run_gapfill = True
        ggf = GeneGapFiller()
        ggf.find_missing_genes(config['gapfilling']['GeneGapFiller parameters']['gff'],
                               current_libmodel)
        ggf.find_missing_reactions(current_model,
                               mail = config['tech-resources']['email'],
                               check_NCBI = config['gapfilling']['GeneGapFiller parameters']['check-NCBI'],
                               fasta = config['general']['protein_fasta'],
                               dmnd_db = config['gapfilling']['GeneGapFiller parameters']['swissprot-dmnd'],
                               swissprot_map = config['gapfilling']['GeneGapFiller parameters']['swissprot-mapping'],
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
        if config['general']['save_all_models']:
            write_model_to_file(current_libmodel, str(Path(dir,'cmpb_out','models',f'{current_libmodel.getId()}_after_Gene_gapfill.xml')))     
            current_modelpath = Path(dir,'cmpb_out','models',f'{current_libmodel.getId()}_after_Gene_gapfill.xml')
            current_model = load_model(current_modelpath,'cobra')
        else:
            write_model_to_file(current_libmodel, str(only_modelpath))
            current_model = load_model(current_modelpath,'cobra')

    # testing
    if run_gapfill:
        between_growth_test(current_model,config,step='after_gapfill')
        between_analysis(current_model, config, step='after_gapfill')

    # ModelPolisher
    ###############
    # @TODO
    # future update
    # currently being revamped 
    # and python access is coming soon

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
        if config['general']['save_all_models']:
            write_model_to_file(current_libmodel, str(Path(dir,'cmpb_out','models',f'{current_libmodel.getId()}_added_KeggPathwayGroups.xml')))
            current_modelpath = Path(dir,'cmpb_out','models',f'{current_libmodel.getId()}_added_KeggPathwayGroups.xml')
        else:
            write_model_to_file(current_libmodel, str(current_modelpath))

    # SBOannotator
    # ------------
    current_libmodel = load_model(str(current_modelpath),'libsbml')
    
    if config['general']['save_all_models']:
        current_libmodel = run_SBOannotator(current_libmodel)
        write_model_to_file(current_libmodel, str(Path(dir,'cmpb_out','models', f'{current_libmodel.getId()}_SBOannotated.xml')))
        current_modelpath = Path(dir,'cmpb_out','models', f'{current_libmodel.getId()}_SBOannotated.xml')
    else:
        current_libmodel = run_SBOannotator(current_libmodel)
        write_model_to_file(current_libmodel, str(current_modelpath))

    current_model = load_model(str(current_modelpath),'cobra')
    between_analysis(current_model,config,step='after_annotation')
    

    # model cleanup
    ###############
    current_model = load_model(str(current_modelpath), 'cobra')
    
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
    if config['general']['save_all_models']:
        write_model_to_file(current_model, str(Path(dir,'cmpb_out','models',f'{current_model.id}_after_duplicate_removal.xml')))
        current_modelpath = Path(dir,'cmpb_out','models',f'{current_model.id}_after_duplicate_removal.xml')
    else:
        write_model_to_file(current_model, str(only_modelpath))
        current_modelpath = only_modelpath
    
    # in-between testing
    between_growth_test(current_model,config,step='after_duplicate_removal')
    between_analysis(current_model,config,step='after_duplicate_removal')

    # BOF
    # ---
    # BOFdat - optional
    current_model = load_model(str(current_modelpath),'cobra')
    if config['BOF']['run_bofdat']:
        check_bof = test_biomass_presence(current_model)
        if check_bof:
            current_model.reactions.get_by_id(check_bof[0]).reaction = adjust_BOF(config['BOF']['bofdat_params']['full_genome_sequence'], 
                                current_modelpath,
                                current_model, 
                                dna_weight_fraction = config['BOF']['bofdat_params']['dna_weight_fraction'], 
                                weight_frac = config['BOF']['bofdat_params']['weight_fraction']) 
            current_model = check_normalise_biomass(current_model)
        else: 
            bof_reac = Reaction('Biomass_BOFdat')
            bof_reac.name = 'Biomass objective function created by BOFdat'
            current_model.add_reactions([bof_reac])
            current_model.reactions.get_by_id(bof_reac).reaction = adjust_BOF(config['BOF']['bofdat_params']['full_genome_sequence'], 
                                current_modelpath,
                                current_model, 
                                dna_weight_fraction = config['BOF']['bofdat_params']['dna_weight_fraction'], 
                                weight_frac = config['BOF']['bofdat_params']['weight_fraction'])
            current_model = check_normalise_biomass(current_model)
    # check and normalise 
    else:
        current_model = check_normalise_biomass(current_model)

    # save model
    if config['general']['save_all_models']:
        write_model_to_file(current_model, str(Path(dir,'cmpb_out','models',f'{current_model.id}_after_BOF.xml')))
        current_modelpath = Path(dir,'cmpb_out','models',f'{current_model.id}_after_BOF.xml')
    else:
        write_model_to_file(current_model, str(current_modelpath))
    
    # MCC
    # ---
    Path(dir,"cmpb_out",'misc','mcc').mkdir(parents=True, exist_ok=False)
    current_model = perform_mcc(current_model, Path(dir,'cmpb_out','misc','mcc'),apply=True)

    # save the final model
    if config['general']['save_all_models']:
        write_model_to_file(current_model, str(Path(dir,'cmpb_out','models',f'{current_model.id}_final_model.xml')))
        current_modelpath = Path(dir,'cmpb_out','models',f'{current_model.id}_final_model.xml')
    else:
        write_model_to_file(current_model, str(only_modelpath))
        current_modelpath = only_modelpath

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
    """Run given settings for the CarveMe-ModelPolisher-based (CMPB) workflow on multuple models.
    """
    pass