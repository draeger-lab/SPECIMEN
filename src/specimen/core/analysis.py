"""Analyse a model (part 5 of the workflow).
"""

__author__ = 'Carolin Brune'

################################################################################
# requirements
################################################################################

import refinegems as rg
import time
import warnings

from pathlib import Path

from .. import util

################################################################################
# functions
################################################################################

# run this part
# -------------
# @NOTE: The below completely changes the parameters, change and check the other connection as needed (e.g. setup, workflow)

def run(model_path, dir, 
        media_path=None, namespace='BiGG',
        pc_model_path=None, pc_based_on='id', 
        test_aa_auxotrophies=True, pathway=True):

    total_time_s = time.time()

    print('\nanalysis\n################################################################################\n')

    # -----------------------
    # create output directory
    # -----------------------

    # make sure given directory path ends with '/'
    if not dir.endswith('/'):
        dir = dir + '/'

    try:
        Path(F"{dir}05_analysis/").mkdir(parents=True, exist_ok=False)
        print(F'Creating new directory {F"{dir}05_analysis/"}')
    except FileExistsError:
        print('Given directory already has required structure.')

    # load model
    model = rg.io.load_model_cobra(model_path)

    # ------------------
    # general statistics
    # ------------------

    print('\n# ------------------\n# general statistics\n# ------------------')

    statistics_report = util.cobra_models.generate_statistics(model)
    statistics_report.save(F'{dir}05_analysis/')

    # -----------------
    # pan-core analysis
    # -----------------

    if pc_model_path:
        print('\n# ------------------\n# pan-core analysis\n# ------------------')
        pc_model = rg.io.load_model_cobra(pc_model_path) 
        pan_core_report = rg.core_pan.compare_to_core_pan(model, pc_model, pc_based_on)
        pan_core_report.save(F'{dir}05_analysis/')

    # ----------------
    # pathway analysis
    # ----------------

    if pathway:
        print('\n# -----------------\n# pathway analysis\n# -----------------')
        pathway_report = rg.pathways.kegg_pathway_analysis(model)
        pathway_report.save(F'{dir}05_analysis/')

    # ---------------
    # growth analysis
    # ---------------

    if media_path:
        print('\n# ---------------\n# growth analysis\n# ---------------')

        # try to set objective to growth
        growth_func_list = rg.biomass.test_biomass_presence(model)
        if growth_func_list:
            # independently of how many growth functions are found, the first one will be used
            model.objective = growth_func_list[0]
            # simulate growth on different media
            growth_report = rg.growth.growth_analysis(model, media_path, namespace=namespace, retrieve='report')
            growth_report.save(F'{dir}05_analysis/')

        else:
            warnings.warn('No growth/biomass function detected, growth simulation will be skipped.')

        # test auxotrophies
        if test_aa_auxotrophies:
            media_list = rg.growth.read_media_config(media_path)
            auxo_report = rg.growth.test_auxotrophies(model, media_list[0], media_list[1], namespace)
            auxo_report.save(F'{dir}05_analysis/')

    total_time_e = time.time()
    print(F'total runtime: {total_time_e-total_time_s}')
