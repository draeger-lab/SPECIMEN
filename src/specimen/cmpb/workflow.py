#!/usr/bin/env python

__author__ = "Tobias Fehrenbach, Famke Baeuerle, Gwendolyn O. Döbel and Carolin Brune"

################################################################################
# requirements
################################################################################

import os
import cobra
import logging
import pandas as pd
from datetime import date
from pathlib import Path

import warnings
import yaml

from libsbml import readSBML

from refinegems.curation.pathways import kegg_pathway_analysis
from refinegems.classes.reports import ModelInfoReport
from refinegems.curation.biomass import test_biomass_presence
from refinegems.analysis import growth
from refinegems.utility.connections import run_memote, perform_mcc, adjust_BOF
from refinegems.curation.curate import resolve_duplicates
from refinegems.curation.pathways import kegg_pathways
from refinegems.utility.io import load_model, write_model_to_file
from refinegems.curation.polish import polish
from refinegems.curation.charges import correct_charges_modelseed

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
        only_modelpath = Path(dir,'cmpb_out','model.xml')

    # create directory structure
    # --------------------------

    dir = config['general']['dir']
    Path(dir,"cmpb_out").mkdir(parents=True, exist_ok=False)                   # cmpb_out
    Path(dir,"cmpb_out",'models').mkdir(parents=True, exist_ok=False)          #   |- models
    Path(dir,"cmpb_out",'logs').mkdir(parents=True, exist_ok=False)            #   |- logs
    Path(dir,"cmpb_out",'misc').mkdir(parents=True, exist_ok=False)            #   |- misc
    Path(dir,"cmpb_out",'misc', 'memote').mkdir(parents=True, exist_ok=False)  #      |- memote
    Path(dir,"cmpb_out",'misc', 'growth').mkdir(parents=True, exist_ok=False)  #      |- growth
    Path(dir,"cmpb_out",'misc', 'stats').mkdir(parents=True, exist_ok=False)   #      |- stats
    
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
        raise ValueError('Currently, CarveMe has not been included in the pipeline. Please use it separatly.mThis wfunction will be provided in a future update.')
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
            write_model_to_file(Path(dir,'cmpb_out','models','after_CarveMe_correction.xml'))
            current_modelpath = Path(dir,'cmpb_out','models','after_CarveMe_correction.xml')
        else:
            write_model_to_file(only_modelpath)
            current_modelpath = only_modelpath


    # growth test
    # -----------
    current_model = load_model(current_modelpath,'cobra')
    between_growth_test(current_model,config,step='after_CarveMe_correction')
    between_analysis(current_model,config,step='after_CarveMe_correction')


    # gapfilling
    ############
    # @TODO
    # options: automatic/manual extension/manual input

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
            write_model_to_file(Path(dir,'cmpb_out','models','added_KeggPathwayGroups.xml'))
            current_modelpath = Path(dir,'cmpb_out','models','added_KeggPathwayGroups.xml')
        else:
            write_model_to_file(current_modelpath)

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
        write_model_to_file(Path(dir,'cmpb_out','models','after_duplicate_removal.xml'))
        current_modelpath = Path(dir,'cmpb_out','models','after_duplicate_removal.xml')
    else:
        write_model_to_file(only_modelpath)
        current_modelpath = only_modelpath

    # BOF
    # ---
    # @TODO
    # BOFdat - optional
    if config['BOF']['run_bofdat']:
        new_bof = adjust_BOF(config['BOF']['bofdat_params']['full_genome_sequence'], 
                            current_modelpath,
                            current_model, 
                            dna_weight_fraction = config['BOF']['bofdat_params']['dna_weight_fraction'], 
                            weight_frac = config['BOF']['bofdat_params']['weight_fraction']) 
        # @TODO 
        # save ?
    # check and normalise
    elif config['BOF']['normalise']:
        pass
        #@TODO
    
    # MCC
    # ---
    # @TODO 
    current_model = perform_mcc(current_model, Path(dir,'cmpb','misc','mcc'),apply=True)

    # save the final model
    if config['general']['save_all_models']:
        write_model_to_file(Path(dir,'cmpb_out','models','final_model.xml'))
        current_modelpath = Path(dir,'cmpb_out','models','final_model.xml')
    else:
        write_model_to_file(only_modelpath)
        current_modelpath = only_modelpath

    # analysis
    ##########
    # @TODO 
    current_model = load_model(current_modelpath,'cobra')
    current_libmodel = load_model(current_modelpath,'libsbml')

    # stats
    # -----
    stats_report = ModelInfoReport(current_model)
    stats_report.save(Path(dir,'stats')) # @TODO adjust Path, just a placeholder really
    
    # kegg pathway
    # ------------
    pathway_report = kegg_pathway_analysis(current_model)
    pathway_report.save(Path(dir,'kegg_pathway')) # @TODO adjust Path, just a placeholder really
    
    # sbo term
    # --------
    # @TODO
    fig = plot_rea_sbo_single(current_libmodel)

    # memote
    # ------
    run_memote(current_model, 'html', save_res=Path(dir,'memote','final_memote.html'))

    # growth
    # ------
    # try to set objective to growth
    growth_func_list = test_biomass_presence(current_model)
    if growth_func_list:
        # independently of how many growth functions are found, the first one will be used
        current_model.objective = growth_func_list[0]
        # simulate growth on different media
        growth_report = growth.growth_analysis(current_model, config['input']['mediapath'], 
                                               namespace=config['general']['namespace'], 
                                               retrieve='report')
        growth_report.save(Path(dir,'growth')) # @TODO adjust Path, just a placeholder really

    else:
        warnings.warn('No growth/biomass function detected, final growth simulation will be skipped.')

    # auxotrophies
    # ------------
    media_list = growth.read_media_config(config['input']['mediapath'])
    auxo_report = growth.test_auxotrophies(current_model, media_list[0], media_list[1], config['general']['namespace'])
    auxo_report.save(Path(dir,'auxotrophies')) # @TODO adjust Path, just a placeholder really

    






###########
# old stuff
###########
# @TODO:
# despite some changes, this is still mainly the old version, 
# which WILL NOT work properly with the new refinegems update - yet
# address TODOs before release!!!!!

def run_old(configpath=None):
    """main function to run the program"""
    print("RefineGEMs provides functions to curate and investigate genome-scale metabolic models!")
    print("Author:", __author__)
    
    config = save_cmpb_user_input(configpath)
    today = date.today().strftime("%Y%m%d")
    log_file = Path(config["out_path"],f'rg_{str(today)}.log')
    
    # check if the output directory is already present, if not create it
    if not os.path.isdir(config['out_path']):
        logging.info('Given out_path is not yet a directory, creating ' + config['out_path'])
        os.makedirs(config['out_path'])
    if config['visualize']:
        dir = Path(config['out_path'],'visualization')
        if not os.path.isdir(dir): 
            os.makedirs(dir) 
    
    print(f'The following logs are saved to: {log_file}')

    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s:%(name)s] %(message)s",
        handlers=[
            logging.FileHandler(log_file),
            logging.StreamHandler()
        ]
    )
    logging.getLogger('cobra').setLevel(logging.WARNING)
    logging.getLogger('requests').setLevel(logging.WARNING)
    logging.getLogger('urllib3').setLevel(logging.WARNING)
    logging.info('----------- New run of refineGEMs -----------')
    
    logging.info(f'Your output will be saved to: {config["out_path"]}')
    
    if config['multiple']:
        logging.info('Growth simulation for multiple models: ')
        models_cobra = rg.utility.io.load_model(config['multiple_paths'], 'cobra')
        growth_all = rg.analysis.growth.growth_analysis(model_cobra,
                    config['mediapath'],
                    namespace = config['namespace'],
                    retrieve='report')
        growth_all.save(Path(config['out_path']))
        logging.info('Multiple model growth simulation results are saved to ' +  str(Path(config['out_path'], 'GrowthSimReport')))

        ###########

        # visualizations
        # .....................................
        # @TODO: more than half of the functions work differently now
        if config['visualize']:
            models_libsbml = rg.io.load_multiple_models(config['multiple_paths'], 'libsbml')
            ini_plot = rg.comparison.plot_initial_analysis(models_libsbml).get_figure()
            sbo_fig_all = rg.comparison.plot_rea_sbo_multiple(models_libsbml).get_figure()
            venn_reac = rg.comparison.plot_venn(models_cobra, 'reaction', True).get_figure()
            venn_metab = rg.comparison.plot_venn(models_cobra, 'metabolite', True).get_figure()
            heatmap = rg.comparison.plot_heatmap_dt(growth_all[['model', 'medium', 'doubling_time [min]']])
            native_heatmap = rg.comparison.plot_heatmap_native(growth_all)
            # saving them
            sbo_fig_all.savefig(config['out_path'] + 'visualization/' + 'all_ReacPerSBO_' + str(today) + '.png', bbox_inches='tight')
            venn_reac.savefig(config['out_path'] + 'visualization/' + 'all_ReacOverlap_' + str(today) + '.png', bbox_inches='tight')
            venn_metab.savefig(config['out_path'] + 'visualization/' + 'all_MetabOverlap_' + str(today) + '.png', bbox_inches='tight')
            heatmap_dt_prefix = 'heatmap_dt_additives_anaerobic_' if config['anaerobic_growth'] else 'heatmap_dt_additives_'
            heatmap.savefig(config['out_path'] + 'visualization/' + heatmap_dt_prefix + str(today) + '.png')
            native_heatmap_prefix = 'heatmap_native_anaerobic_' if config['anaerobic_growth'] else 'heatmap_native_'
            native_heatmap.savefig(config['out_path'] + 'visualization/' + native_heatmap_prefix + str(today) + '.png', bbox_inches='tight')
            ini_plot.savefig(config['out_path'] + 'visualization/' + 'model_status_' + str(today) + '.png', bbox_inches='tight')
        # .....................................

    if config['single']:        
        try:    
            model_cobra, errors = cobra.io.sbml.validate_sbml_model(config['model'])
            #logging.info(errors)
        except (OSError):
            model_cobra = None
            logging.info('No or no valid model given, please enter a valid path in the model field in the config file.')

        if config['keggpathways']:
            model_libsbml, non_kegg = rg.curation.pathways.kegg_pathways(config['model'])
            file = open(config['out_path'] + model_libsbml.getId() + '_reac_wo_kegg_' + str(today) + '.txt','w')
            for reaction in non_kegg:
                file.write(reaction+"\n")
            file.close()
            logging.info('Kegg Pathways were added to the model as groups. Reactions that have no KEGG annotation are denoted in ' + model_libsbml.getId() + '_reac_wo_kegg.txt')

        else:
            model_libsbml = rg.utility.io.load_model(config['model'],'libsbml')
        
        if config['sboterms']:
            if config['visualize']:
                # .....................................
                # @TODO: where us the functions now?
                sbo_fig = rg.investigate.plot_rea_sbo_single(model_libsbml).get_figure()
                # .....................................
                # saving the created visualizations
                sbo_fig.savefig(config['out_path'] + 'visualization/' + str(model_cobra.id) + '_ReacPerSBO_beforeUpdate_' + str(today) + '.png', bbox_inches='tight')
            model_libsbml = rg.sboann.sbo_annotation(model_libsbml)
            logging.info('SBO Terms updated for ' + model_libsbml.getId())
            
        if config['charge_corr']:
            model_libsbml, multiple_charges = rg.curation.charges.correct_charges_modelseed(model_libsbml)
            pd.DataFrame.from_dict(multiple_charges, orient='index').to_csv(config['out_path'] + model_libsbml.getId() + '_mulchar_' + str(today) + '.csv', sep=',', header=False)
            logging.info('Charges were corrected for ' + model_libsbml.getId() + '. A table with metabolites with multiple charges can be found under ' + model_libsbml.getId() + '_mulchar_' + str(today) + '.csv')
            
        if config['man_cur']:
            if config['man_cur_type'] == 'gapfill':
                gapfill = rg.utility.io.load_manual_gapfill(config['man_cur_table'])
                model_libsbml = rg.curation.curate.add_reactions_from_table(model_libsbml, gapfill, config['entrez_email'])
                logging.info('Manual gap filling was done for ' + model_libsbml.getId())
            elif config['man_cur_type'] == 'metabs':
                man_ann = rg.utility.io.load_manual_annotations(config['man_cur_table'])
                model_libsbml = rg.curation.curate.update_annotations_from_table(model_libsbml, man_ann)
                model_libsbml = rg.curation.curate.update_annotations_from_others(model_libsbml)
                logging.info('Manual update of annotations was done for ' + model_libsbml.getId())
                
        if config['gap_analysis'] and config['gapfill_model']:
            filename = f'{config["out_path"]}{model_libsbml.getId()}_gap_analysis_{str(today)}'
            if config['gap_analysis_params'].get('db_to_compare') not in ('BioCyc', 'KEGG+BioCyc'):
                logging.warning('Currently, only the result from the \'BioCyc\' or \'KEGG+BioCyc\' runs can be directly added to a model.')
                gap_analysis = rg.curation.gapfill.gap_analysis(model_libsbml, config['gap_analysis_params'], filename)
                logging.info(f'Gap analysis for {model_libsbml.getId()} with {config["gap_analysis_params"].get("db_to_compare")} was performed.')
                logging.info(f'Complete Excel table is in file: {filename}.')
            else:
                gapfill = rg.curation.gapfill.gapfill(model_libsbml, config['gap_analysis_params'], filename)
                gap_analysis_stats = gapfill[0][0]
                logging.info(f'Statistics on missing entites for {model_libsbml.getId()}:')
                logging.info(gap_analysis_stats)
                logging.info(f'Complete Excel table is in file: {filename}.')
                model_libsbml = gapfill[-1]
                logging.info(f'Gaps were filled in {model_libsbml.getId()}.')
        elif config['gap_analysis']:
            filename = f'{config["out_path"]}{model_libsbml.getId()}_gap_analysis_{str(today)}'
            gap_analysis = rg.curation.gapfill.gap_analysis(model_libsbml, config['gap_analysis_params'], filename)
            logging.info(f'Gap analysis for {model_libsbml.getId()} with {config["gap_analysis_params"].get("db_to_compare")} was performed.')
            if  config["gap_analysis_params"].get("db_to_compare") != 'KEGG':
                logging.info(f'Statistics on missing entites for {model_libsbml.getId()}:')
                logging.info(gap_analysis[0])
            logging.info(f'Complete Excel table is in file: {filename}.')
        elif config['gapfill_model']:
            model_libsbml = rg.curation.gapfill.gapfill_model(model_libsbml, config['gap_analysis_file'])
            logging.info(f'Gaps were filled in {model_libsbml.getId()}.')
 
        if config['polish']:
            model_libsbml = rg.curation.polish.polish(model_libsbml, config['entrez_email'], config['id_db'], config['gff_file'], 
                                             config['protein_fasta'], config['lab_strain'], config['organismid'], config['out_path'])
            logging.info(model_libsbml.getId() + ' has been polished')
            
        if config['biomass']:
            result = rg.curation.biomass.check_normalise_biomass(model_cobra)
            if result:
                model_libsbml = result
                logging.info(model_libsbml.getId() + '\'s biomass function has been checked.')
        
        mods = [config['keggpathways'], config['sboterms'], config['charge_corr'], config['man_cur'], config['gapfill_model'], config['polish'], config['biomass']]
        
        if any(mods):
            if config['model_out'] == 'stdout':   
                config['model_out'] = config['out_path'] + model_libsbml.getId() + '_modified_' + str(today) + '.xml'

            try:
                rg.utility.io.write_model_to_file(model_libsbml, config['model_out'])
            except (OSError) as e:
                logging.info(e)
                logging.info("Model could not be saved...")
            
            if model_cobra:                                          
                try:    
                    model_cobra, errors = cobra.io.sbml.validate_sbml_model(config['model_out'])
                    logging.info(errors)
                except (OSError) as e:
                    logging.info(e)
                    model_cobra = None
                    logging.info('Model was invalidated during curation steps.')

        if model_cobra:
            logging.info(model_cobra.id + ' will be investigated.')
            # .....................................................
            # @TODO: different !!!!!
            name, reac, metab, genes = rg.analysis.investigate.initial_analysis(model_libsbml) 
            orphans, deadends, disconnected = rg.investigate.get_orphans_deadends_disconnected(model_cobra)
            mass_unbal, charge_unbal = rg.investigate.get_mass_charge_unbalanced(model_cobra)
            egc = rg.investigate.get_egc(model_cobra)
            if config['visualize']:
                logging.info('All visualizations can be found in the subfolder "visualization".')
                sbo_fig = rg.investigate.plot_rea_sbo_single(model_libsbml).get_figure()
                
                # saving the created visualizations
                sbo_fig.savefig(config['out_path'] + 'visualization/' + str(model_cobra.id) + '_ReacPerSBO_' + str(today) + '.png', bbox_inches='tight')
            # .....................................................
            # .....................................................
            # @TODO different function
            if config['memote']:
                score = rg.investigate.get_memote_score(rg.investigate.run_memote(model_cobra))
            # .....................................................
            if config['modelseed']:
                charge_mismatch, formula_mismatch = rg.curation.db_access.modelseed.compare_to_modelseed(model_cobra)
            # .....................................................
            # @TODO function changed
            if config['media']:
                growth_sim = rg.growth.get_growth_selected_media(model_cobra, config['media'], config['growth_basis'], config['anaerobic_growth'])
            # .....................................................
            if config['memote']:
                information = [[name], [reac], [metab], [genes], [score], orphans, deadends, disconnected, mass_unbal, charge_unbal]
                model_params = pd.DataFrame(information, ['model name', '#reactions', '#metabolites', '#genes', 'memote score', 'orphans', 'deadends', 'disconnected', 'mass unbalanced', 'charge unbalanced']).T
            else:
                information = [[name], [reac], [metab], [genes], orphans, deadends, disconnected, mass_unbal, charge_unbal]
                model_params = pd.DataFrame(information, ['model name', '#reactions', '#metabolites', '#genes', 'orphans', 'deadends', 'disconnected', 'mass unbalanced', 'charge unbalanced']).T
            with pd.ExcelWriter(config['out_path'] + name + '_' + str(today) + '.xlsx') as writer:  
                model_params.to_excel(writer, sheet_name='model params', index=False)
                growth_sim_name = 'anaerobic growth simulation' if config['anaerobic_growth'] else 'growth simulation'
                if config['media']: #@TODO: needed???
                    growth_sim.to_excel(writer, sheet_name=growth_sim_name, index=False)
                egc.to_excel(writer, sheet_name='EGC test', index=False)
                if config['modelseed']:
                    charge_mismatch.to_excel(writer, sheet_name='charge mismatches', index=False)
                    formula_mismatch.to_excel(writer, sheet_name='formula mismatches', index=False)
                logging.info('Single model analyses results are saved to ' + name + '_analyses_' + str(today) + '.xlsx')
        else:
            logging.info('No valid model, investigation aborted!')


####### IDEAS below ####

# run for multiple models
def wrapper():
    pass