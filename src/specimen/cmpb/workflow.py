#!/usr/bin/env python

__author__ = "Tobias Fehrenbach, Famke Baeuerle, Gwendolyn O. Döbel and Carolin Brune"

################################################################################
# requirements
################################################################################

import logging
import pandas as pd
from datetime import date
from pathlib import Path

import warnings
import yaml

from cobra import Reaction
from libsbml import readSBML

from refinegems.analysis import growth
from refinegems.analysis.investigate import plot_rea_sbo_single
from refinegems.classes.reports import ModelInfoReport
from refinegems.curation.biomass import test_biomass_presence, check_normalise_biomass
from refinegems.curation.charges import correct_charges_modelseed
from refinegems.curation.curate import resolve_duplicates
from refinegems.curation.gapfill import gapfill_model, gapfill
from refinegems.curation.pathways import kegg_pathways, kegg_pathway_analysis
from refinegems.curation.polish import polish
from refinegems.utility.connections import run_memote, perform_mcc, adjust_BOF
from refinegems.utility.io import load_model, write_model_to_file

# from SBOannotator import *
from SBOannotator import sbo_annotator

from ..util.set_up import save_cmpb_user_input

################################################################################
# functions
################################################################################

    # ....................................................
    # @TODO / @IDEAS
    # global options
    # run memote after every step
    # calculate model stats after each step
    # use temp folder or report all model/in-between steps
    # what to write in the log file
    # ....................................................

# dev notes
#   in the run function: current_model means the cobrapy model, 
#   while current_libmodel means the libsbml model

# @TODO / @IDEAS Add option to have specific colour list per model for plots
# @TODO Maybe get models at first and then add model IDs to every save filename?

def run(configpath:str):

    def between_growth_test(model, cfg, step):
        # try to set objective to growth
        growth_func_list = test_biomass_presence(model)
        if growth_func_list:
            # independently of how many growth functions are found, the first one will be used
            model.objective = growth_func_list[0]
            # simulate growth on different media
            growth_report = growth.growth_analysis(model, cfg['input']['mediapath'], 
                                                namespace=cfg['general']['namespace'], 
                                                retrieve='report')
            growth_report.save(Path(cfg['general']["dir"], 'cmpb_out', 'misc', 'growth', step)) 
        else:
            mes = f'No growth/biomass function detected, growth simulation for step {step} will be skipped.'
            warnings.warn(mes)

    def between_analysis(model, cfg, step):
        # optional analysis
        if cfg['general']['memote_always_on']:
            run_memote(model, 'html', 
                    return_res=False, 
                    save_res=Path(cfg['general']["dir"], 'cmpb_out', 'misc','memote',f'{step}_.html'), 
                    verbose=False)
        if cfg['general']['stats_always_on']:
            report = ModelInfoReport(model)
            report.save(Path(cfg['general']["dir"],'cmpb_out', 'misc', 'stats',f'{step}_.html')) 


    # setup phase
    #############

    # load config
    # -----------
    if not configpath:
        config = save_cmpb_user_input(Path('cmpb_out', 'logs', 'config_user.yaml')) 
    else:
        with open(configpath, "r") as cfg:
            config = yaml.load(cfg, Loader=yaml.loader.FullLoader)

    if not config['general']['save_all_models']:
        only_modelpath = Path(dir,'cmpb_out','model.xml') # @TODO Use model ID here...

    # create directory structure
    # --------------------------

    dir = config['general']['dir']
    Path(dir,"cmpb_out").mkdir(parents=True, exist_ok=False)                          # cmpb_out
    Path(dir,"cmpb_out",'models').mkdir(parents=True, exist_ok=False)                 #   |- models
    Path(dir,"cmpb_out",'logs').mkdir(parents=True, exist_ok=False)                   #   |- logs
    Path(dir,"cmpb_out",'misc').mkdir(parents=True, exist_ok=False)                   #   |- misc
    Path(dir,"cmpb_out",'misc', 'memote').mkdir(parents=True, exist_ok=False)         #      |- memote
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
    # @TODO
    # will come in a future update
    if not config['input']['modelpath']:
        # run CarveMe
        raise ValueError('Currently, CarveMe has not been included in the pipeline. Please use it separatly. This function will be provided in a future update.')
    else:
        current_modelpath = config['input']['modelpath']

    # optional analysis
    current_model = load_model(current_modelpath,'cobra')
    between_analysis(current_model, config, step='after_draft')

    # CarveMe correction
    ####################

    current_libmodel = load_model(current_modelpath,'libsbml')
    # check, if input is a CarveMe model
    if 'CarveMe' in current_libmodel.getNotesString():
        Path(dir,"cmpb_out",'misc', 'wrong_annotations').mkdir(parents=True, exist_ok=False)
        current_libmodel = polish(current_libmodel, 
                                  email = config['cm-polish']['email'], 
                                  id_db = config['general']['namespace'], 
                                  refseq_gff = config['cm-polish']['refseq_gff'], 
                                  protein_fasta = config['cm-polish']['protein_fasta'], 
                                  lab_strain = config['cm-polish']['is_lab_strain'], 
                                  kegg_organism_id = config['cm-polish']['kegg_organism_id'], 
                                  path = Path(dir,'cmpb_out','misc','wrong_annotations')) 
        # rg correct charges
        current_libmodel, mult_charges_dict = correct_charges_modelseed(current_libmodel)
        mult_charges_tab = pd.DataFrame.from_dict(mult_charges_dict, orient='index')
        mult_charges_tab.to_csv(Path(dir,'cmpb_out','misc','reac_with_mult_charges.tsv'), sep='\t')
        
        # save model
        if config['general']['save_all_models']:
            write_model_to_file(current_libmodel, Path(dir,'cmpb_out','models',f'{current_libmodel.getId()}_after_CarveMe_correction.xml'))
            current_modelpath = Path(dir,'cmpb_out','models',f'{current_libmodel.getId()}_after_CarveMe_correction.xml')
        else:
            write_model_to_file(current_libmodel, only_modelpath)
            current_modelpath = only_modelpath


    # growth test
    # -----------
    current_model = load_model(current_modelpath,'cobra')
    between_growth_test(current_model,config,step='after_CarveMe_correction')
    between_analysis(current_model,config,step='after_CarveMe_correction')


    # gapfilling
    ############
    # options: automatic/manual extension/manual input
    if config['gapfilling']['gap_fill_params']: 
        
        filename = Path(dir, 'cmpb', 'misc',f'{current_libmodel.getId()}_gap_analysis_{str(today)}')
        
        match config['gapfilling']['gap_fill_params']['db_to_compare']:
            
            case 'KEGG':
                gap_analysis_result, current_libmodel = gapfill(
                    current_libmodel, 'KEGG', config['general']['refseq_gff'], 
                    config['general']['kegg_organism_id'], None, 
                    filename
                    )
            case 'BioCyc':
                gap_analysis_result, current_libmodel = gapfill(
                    current_libmodel, 'BioCyc', None, None, config['gapfilling']['gap_fill_params']['biocyc_files'], 
                    filename
                    )
            case 'KEGG+BioCyc':
                gap_analysis_result, current_libmodel = gapfill(
                    current_libmodel, 'KEGG+BioCyc', config['general']['refseq_gff'],
                    config['general']['kegg_organism_id'], config['gapfilling']['gap_fill_params']['biocyc_files'], 
                    filename
                    )
            case _:
                logging.warning(
                    f'''
                    \'{config["gapfilling"]["gap_fill_params"]["db_to_compare"]}\' is no valid input for the gap analysis.
                    Please specify for the parameter \'db_to_compare\' one of the following options: 
                    KEGG | BioCyc | KEGG+BioCyc
                    '''
                    )
                
        logging.info(f'Gap analysis for {current_libmodel.getId()} with {config['gapfilling']['gap_fill_params']['db_to_compare']} was performed.')
        logging.info(f'Complete Excel table is in file: {filename}.xlsx.')

        # save model
        if config['general']['save_all_models']:
            write_model_to_file(current_libmodel, Path(dir,'cmpb_out','models',f'{current_libmodel.getId()}_autofilled_gaps.xml'))
            current_modelpath = Path(dir,'cmpb_out','models',f'{current_libmodel.getId()}_autofilled_gaps.xml')
        else:
            write_model_to_file(current_libmodel, current_modelpath)
            
        logging.info(f'Gaps were filled automatically in {current_libmodel.getId()}.')
                
    if config['gapfilling']['gap_fill_file']: 
        current_libmodel = gapfill_model(current_libmodel, config['gapfilling']['gap_fill_file'])

        # save model
        if config['general']['save_all_models']:
            write_model_to_file(current_libmodel, Path(dir,'cmpb_out','models',f'{current_libmodel.getId()}_manually_filled_gaps.xml'))
            current_modelpath = Path(dir,'cmpb_out','models',f'{current_libmodel.getId()}_manually_filled_gaps.xml')
        else:
            write_model_to_file(current_libmodel, current_modelpath)
            
        logging.info(f'Gaps were filled based on a  manually curated file in {current_libmodel.getId()}.')
        

    # ModelPolisher
    ###############
    # @TODO

    # Annotations
    #############

    # KEGGPathwayGroups
    # -----------------
    if config['kegg_pathway_groups']:
        current_libmodel, missing_list = kegg_pathways(current_modelpath)
        with open(Path(dir, 'cmpb', 'misc', 'reac_wo_kegg_pathway_groups.txt'), 'w') as outfile:
            for line in missing_list:
                outfile.write(f'{line}\n')
        # save model
        if config['general']['save_all_models']:
            write_model_to_file(current_libmodel, Path(dir,'cmpb_out','models',f'{current_libmodel.getId()}_added_KeggPathwayGroups.xml'))
            current_modelpath = Path(dir,'cmpb_out','models',f'{current_libmodel.getId()}_added_KeggPathwayGroups.xml')
        else:
            write_model_to_file(current_libmodel, current_modelpath)

    # SBOannotator
    # ------------
    # @TODO 
    #  theoretically: something along the way:
    libsbml_doc = readSBML(current_modelpath)
    libsbml_model = libsbml_doc.getModel()
    if config['general']['save_all_models']:
        sbo_annotator(libsbml_doc, libsbml_model, 'constraint-based', True, 'create_dbs', 
                      Path(dir,'cmpb','models', 'SBOannotated.xml'))
    else:
        sbo_annotator(libsbml_doc, libsbml_model, 'constraint-based', True, 'create_dbs', 
                      Path(current_modelpath))

    between_analysis(current_model,config,step='after_annotation')
    

    # model cleanup
    ###############
    current_model = load_model(current_modelpath)

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

    between_growth_test(current_model,config,step='after_duplicate_removal')
    between_analysis(current_model,config,step='after_duplicate_removal')

    # save model
    if config['general']['save_all_models']:
        write_model_to_file(current_model, Path(dir,'cmpb_out','models',f'{current_model.getId()}_after_duplicate_removal.xml'))
        current_modelpath = Path(dir,'cmpb_out','models',f'{current_model.getId()}_after_duplicate_removal.xml')
    else:
        write_model_to_file(current_model, only_modelpath)
        current_modelpath = only_modelpath

    # BOF
    # ---
    # BOFdat - optional
    current_model = load_model(current_modelpath,'cobra')
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

    # save 
    # save model
    if config['general']['save_all_models']:
        write_model_to_file(current_model, Path(dir,'cmpb_out','models',f'{current_model.getId()}_after_BOF.xml'))
        current_modelpath = Path(dir,'cmpb_out','models',f'{current_model.getId()}_after_BOF.xml')
    else:
        write_model_to_file(current_model, current_modelpath)
    
    # MCC
    # ---
    current_model = perform_mcc(current_model, Path(dir,'cmpb','misc','mcc'),apply=True)

    # save the final model
    if config['general']['save_all_models']:
        write_model_to_file(current_model, Path(dir,'cmpb_out','models',f'{current_model.getId()}_final_model.xml'))
        current_modelpath = Path(dir,'cmpb_out','models',f'{current_model.getId()}_final_model.xml')
    else:
        write_model_to_file(current_model, only_modelpath)
        current_modelpath = only_modelpath

    # analysis
    ##########
    current_model = load_model(current_modelpath,'cobra')
    current_libmodel = load_model(current_modelpath,'libsbml')

    # stats
    # -----
    stats_report = ModelInfoReport(current_model)
    stats_report.save(Path(dir,'cmpb_out', 'misc','stats')) 
    
    # kegg pathway
    # ------------
    pathway_report = kegg_pathway_analysis(current_model)
    pathway_report.save(Path(dir,'cmpb_out','misc','kegg_pathway')) 
    
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
    auxo_report.save(Path(dir,'cmpb_out','misc','auxotrophies')) 


####### IDEAS below ####

# run for multiple models
def wrapper():
    pass