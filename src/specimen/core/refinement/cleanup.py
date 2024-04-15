"""Perform part 3 of the refinement: clean-up.
"""

__author__ = 'Carolin Brune'

################################################################################
# requirements
################################################################################

import cobra
import math
import numpy as np
import pandas as pd
import sys
import time
import warnings

from pathlib import Path
from typing import Literal,Union

from refinegems.utility.io import load_model
from refinegems.classes.medium import medium_to_model, Medium
from refinegems.analysis.growth import read_media_config
from refinegems.analysis.investigate import run_memote
from refinegems.curation.curate import resolve_duplicates, complete_BioMetaCyc

################################################################################
# variables
################################################################################

################################################################################
# functions
################################################################################


# @ADDED TO refinegems 
# @TODO - use the functions here as well, without skipping the creation label
# --> how to include the check? 
# - for now, this func will be kept
def check_direction(model:cobra.Model,data_file:str) -> cobra.Model:
    """Check the direction of newly created reactions (01-extention) by searching for matching MetaCyc,
    KEGG and MetaNetX IDs as well as EC number in a downloaded BioCyc (MetaCyc)
    database table (need to contain at least the following columns:
    Reactions (MetaCyc ID),EC-Number,KEGG reaction,METANETX,Reaction-Direction.

    @NOTE: checks only creations that do not contain the notes['creation'] == 'via template'

    Args:
        - model (cobra.Model): The GEM containing the reactions to be check for direction.
        - data_file (str): Path to the MetabCyc (BioCyc) smart table.

    Returns:
        cobra.Model: The updated model.
    """

    # create MetaCyc table
    # --------------------
    data = pd.read_csv(data_file, sep='\t')
    # rewrite the columns into a better comparable/searchable format
    data['KEGG reaction'] = data['KEGG reaction'].str.extract('.*>(R\d*)<.*')
    data['METANETX']      = data['METANETX'].str.extract('.*>(MNXR\d*)<.*')
    data['EC-Number']     = data['EC-Number'].str.extract('EC-(.*)')

    # check direction
    # --------------------
    for r in model.reactions:
            # entry from template, assumed to be already curated
            if 'creation' in r.notes and 'via template' == r.notes['creation']:
                continue
            # newly created entry, check direction with BioCyc
            else:
                direction = None
                # easy case: metacyc is already (corretly) annotated
                if 'metacyc.reaction' in r.annotation and len(data[data['Reactions'] == r.annotation['metacyc.reaction']]) != 0:
                    direction = data[data['Reactions'] == r.annotation['metacyc.reaction']]['Reaction-Direction'].iloc[0]
                    r.notes['BioCyc direction check'] = F'found {direction}'
                # complicated case: no metacyc annotation
                else:
                    annotations = []

                    # collect matches
                    if 'kegg.reaction' in r.annotation and r.annotation['kegg.reaction'] in data['KEGG reaction'].tolist():
                        annotations.append(data[data['KEGG reaction'] == r.annotation['kegg.reaction']]['Reactions'].tolist())
                    if 'metanetx.reaction' in r.annotation and r.annotation['metanetx.reaction'] in data['METANETX'].tolist():
                        annotations.append(data[data['METANETX'] == r.annotation['metanetx.reaction']]['Reactions'].tolist())
                    if 'ec-code' in r.annotation and r.annotation['ec-code'] in data['EC-Number'].tolist():
                        annotations.append(data[data['EC-Number'] == r.annotation['ec-code']]['Reactions'].tolist())

                    # check results
                    # no matches
                    if len(annotations) == 0:
                        r.notes['BioCyc direction check'] = 'not found'

                    # matches found
                    else:
                        # built intersection
                        intersec = set(annotations[0]).intersection(*annotations)
                        # case 1: exactly one match remains
                        if len(intersec) == 1:
                            entry = intersec.pop()
                            direction = data[data['Reactions'] == entry]['Reaction-Direction'].iloc[0]
                            r.annotation['metacyc.reaction'] = entry
                            r.notes['BioCyc direction check'] = F'found {direction}'

                        # case 2: multiple matches found -> inconclusive
                        else:
                            r.notes['BioCyc direction check'] = F'found, but inconclusive'

                # update direction if possible and needed
                if not pd.isnull(direction):
                    if 'REVERSIBLE' in direction:
                        # set reaction as reversible by setting default values for upper and lower bounds
                        r.lower_bound = -1000.
                    elif 'RIGHT-TO-LEFT' in direction:
                        # invert the default values for the boundaries
                        r.lower_bound = -1000.
                        r.upper_bound = 0.
                    else:
                        # left to right case is the standart for adding reactions
                        # = nothing left to do
                        continue
    return model


# gap filling
# -----------

def single_cobra_gapfill(model:cobra.Model, 
                         universal:cobra.Model, medium:Medium, 
                         namespace:Literal['BiGG']='BiGG', 
                         growth_threshold:float = 0.05) -> Union[list[str],bool]:
    """Attempt gapfilling (with COBRApy) for a given model to allow growth on a given
    medium.

    Args:
        - model (cobra.Model): The model to perform gapfilling on.
        - universal (cobra.Model):  A model with reactions to be potentially used for the gapfilling.
        - medium (Medium): A medium the model should grow on.
        - namespace (Literal['BiGG'], optional): Namespace to use for the model. 
            Defaults to 'BiGG'.
        - growth_threshold (float, optional):  Minimal rate for the model to be considered growing.
            Defaults to 0.05.

    Returns:
        Union[list[str],True]: List of reactions to be added to the model to allow growth
        or True, if the model already grows.
    """
    
    # perform the gapfilling
    solution = []
    with model as model_copy:
        # set medium model should grow on
        medium_to_model(model_copy, medium, namespace, double_o2=False) 
        # if model does not show growth (depending on threshold), perform gapfilling
        if model_copy.optimize().objective_value < growth_threshold:
            try:
                solution = cobra.flux_analysis.gapfill(model_copy, universal,
                                                       lower_bound = growth_threshold,
                                                       demand_reactions=False)
            except cobra.exceptions.Infeasible:
                warnings.warn(F'Gapfilling for medium {medium.name} failed. Manual curation required.')
        else:
            print(F'Model already grows on medium {medium.name} with objective value of {model_copy.optimize().objective_value}')
            return True

    return solution


def cobra_gapfill_wrapper(model:cobra.Model, universal:cobra.Model, 
                          medium:Medium, namespace:Literal['BiGG']='BiGG',
                          growth_threshold:float = 0.05,
                          iterations:int=3, chunk_size:int=10000) -> cobra.Model:
    """Wrapper for single_cobra_gapfill().

    Either use the full set of reactions in universal model by setting iteration to
    0 or None or use them in randomized chunks for faster runtime (especially useful
    on laptops). Note: when using the second option, be aware that this does not
    test all reaction combinations exhaustively (heuristic approach!!!).

    Args:
        - model (cobra.Model): The model to perform gapfilling on.
        - universal (cobra.Model):  A model with reactions to be potentially used for the gapfilling.
        - medium (Medium): A medium the model should grow on.
        - namespace (Literal['BiGG'], optional): Namespace to use for the model. 
            Options include 'BiGG'.
            Defaults to 'BiGG'.
        - growth_threshold (float, optional): Growth threshold for the gapfilling. 
            Defaults to 0.05.
        - iterations (int, optional): Number of iterations for the heuristic version of the gapfilling.
            If 0 or None is given, uses full set of reactions.
            Defaults to 3. 
        - chunk_size (int, optional): Number of reactions to be used for gapfilling at the same time. 
            If None or 0 is given, use full set, not heuristic.
            Defaults to 10000.

    Returns:
        cobra.Model: The gapfilled model, if a solution was found.
    """

    solution = []

    # run a heuristic approach:
    #     for a given number of iterations, use a subset (chunk_size) of
    #     reactions for the gapfilling
    if (iterations and iterations > 0) and (chunk_size and chunk_size > 0):

        num_reac = len(model.reactions)
        # for each iteration
        for i in range(iterations):
            not_used = [_ for _ in range(0,num_reac)]

            # divide reactions in random subsets
            for chunk in range(math.ceil(num_reac/chunk_size)):
                if len(not_used) > chunk_size:
                    rng = np.random.default_rng()
                    reac_set = rng.choice(not_used, size=chunk_size, replace=False, shuffle=False)
                    not_used = [_ for _ in not_used if _ not in reac_set]
                else:
                    reac_set = not_used

                # get subset of reactions
                subset_reac = cobra.Model('subset_reac')
                for n in reac_set:
                    subset_reac.add_reactions([universal.reactions[n].copy()])

                # gapfilling
                solution = single_cobra_gapfill(model, subset_reac, medium, namespace, growth_threshold)

                if (isinstance(solution,bool) and solution) or (isinstance(solution, list) and len(solution) > 0):
                    break

            if (isinstance(solution,bool) and solution) or (isinstance(solution, list) and len(solution) > 0):
                break

    # use the whole reactions content of the universal model at once
    #     not advised for Laptops and PCs with small computing power
    #     may take a long time
    else:
        solution = single_cobra_gapfill(model, universal, medium, namespace, growth_threshold)

    # if solution is found add new reactions to model
    if isinstance(solution, list) and len(solution) > 0:
        for reac in solution[0]:
            reac.notes['creation'] = 'via gapfilling'
        print(F'Adding {len(solution[0])} reactions to model to ensure growth on medium {medium.name}.')
        model.add_reactions(solution[0])

    return model


def multiple_cobra_gapfill(model: cobra.Model, universal:cobra.Model, 
                           media_list:list[Medium], 
                           namespace:Literal['BiGG']='BiGG', 
                           growth_threshold:float = 0.05, iterations:int=3, chunk_size:int=10000) -> cobra.Model:
    """Perform single_cobra_gapfill() on a list of media.

    Args:
        - model (cobra.Model): The model to be gapfilled.
        - universal (cobra.Model): The model to use reactions for gapfilling from.
        - media_list (list[Medium]): List ofmedia the model is supposed to grow on.
        - growth_threshold (float, optional): Growth threshold for the gapfilling. 
            Defaults to 0.05.
        - iterations (int, optional): Number of iterations for the heuristic version of the gapfilling.
            If 0 or None is given, uses full set of reactions.
            Defaults to 3. 
        - chunk_size (int, optional): Number of reactions to be used for gapfilling at the same time. 
            If None or 0 is given, use full set, not heuristic.
            Defaults to 10000.
    Returns:
        cobra.Model: The gapfilled model, if a solution was found.
    """
    

    for medium in media_list:
        model = cobra_gapfill_wrapper(model,universal,medium, namespace, iterations, chunk_size, growth_threshold)

    return model


# run
# ---
# @TODO add gapfilling from refinegems and more
# @TODO maybe a new / different way to save/store/add the params for gapfilling, e.g. a config?
def run(model:str, dir:str, 
        biocyc_db:str=None, 
        check_dupl_reac:bool = False,
        check_dupl_meta:bool = 'default',
        remove_unused_meta:bool = False, 
        remove_dupl_reac:bool = False, 
        remove_dupl_meta:bool = False,
        universal:str = None, 
        media_path:str = None, 
        namespace:Literal['BiGG']='BiGG', 
        growth_threshold:float = 0.05,
        iterations:int=3, chunk_size:int=10000,
        memote:bool = False):
    """Perform the second refinement step, cleanup, on a model.

    The second refinement step resolves the following issues:
    (1) (optional) checking direction of reactions with BioCyc
    (2) resolve BioCyc/MetaCyc annotation inconsistencies
    (3) find and/or resolve duplicates (reactions and metabolites)
    (4) gapfilling using cobra + universal model with reactions + a set of media (@TODO: under developement)

    Args:
        - model (str): The Path to an sbml model.
        - dir (str): Path to the directory of the output.
        - biocyc_db (str, optional): Path to the BioCyc/MetaCyc reaction database file. 
            Defaults to None, which leads to skipping the direction check.
        - check_dupl_reac (bool, optional): Option to check for duplicate reactions. 
            Defaults to False.
        - check_dupl_meta (bool, optional): Option to check for duplicate metabolites. 
            Defaults to 'default', which checks based on MetaNetX first.
            Further options include 'skip' and 'exhaustive' (check for all possibilities).
        - remove_unused_meta (bool, optional): Option to remove unused metabolites from the model. 
            Defaults to False.
        - remove_dupl_reac (bool, optional): Option to remove the duplicate reactions. 
            True is only applied, if check_dupl_reac is also True.
            Defaults to False.
        - remove_dupl_meta (bool, optional): Option to remove the duplicate metabolites. 
            True is only applied, if check_dupl_meta is also True.
            Defaults to False.
        - universal (str, optional): Path to a universal model for gapfilling. 
            Defaults to None, which skips the gapfilling.
        - media_path (str, optional): Path to a medium config file for gapfilling. 
            Defaults to None.
        - namespace (Literal['BiGG'], optional): Namespace to use for the model. 
            Options include 'BiGG'.
            Defaults to 'BiGG'.
        - growth_threshold (float, optional): Growth threshold for the gapfilling. 
            Defaults to 0.05.
        - iterations (int, optional): Number of iterations for the heuristic version of the gapfilling.
            If 0 or None is given, uses full set of reactions.
            Defaults to 3. 
        - chunk_size (int, optional): Number of reactions to be used for gapfilling at the same time. 
            If None or 0 is given, use full set, not heuristic.
            Defaults to 10000.
        - memote (bool, optional): Option to run memote on the cleaned model. 
            Defaults to False.

    Raises:
        ValueError: Unknown option for check_dupl_meta
    """

    # -------------------
    # checking parameters
    # -------------------
    if not check_dupl_meta in ['default','skip','exhaustive']:
        raise ValueError('Unknown option {check_dupl_meta} for checking duplicate metabolite. Use one of: default, skip, exhaustive')

    # -------------
    # start program
    # -------------
    print('\nrefinement step 2: clean-up\n################################################################################\n')

    # -----------------------
    # create output directory
    # -----------------------

    try:
        Path(dir,"step2-clean-up").mkdir(parents=True, exist_ok=False)
        print(F'Creating new directory {str(Path(dir,"step2-clean-up"))}')
    except FileExistsError:
        print('Given directory already has required structure.')

    model = cobra.io.read_sbml_model(model)

    # --------------------
    # check direction
    # --------------------

    if biocyc_db:

        print('\n# --------------------\n# check direction\n# --------------------')

        start = time.time()

        # check direction
        # @TODO : check namespace independency
        model = check_direction(model,biocyc_db)

        end = time.time()
        print(F'\ttime: {end - start}s')

    # ----------------------------------
    # complete BioCyc/MetaCyc annotation
    # ----------------------------------

    print('\n# ----------------------------------\n# complete BioCyc/MetaCyc annotation\n# ----------------------------------')
    start = time.time()

    model = complete_BioMetaCyc(model)

    end = time.time()
    print(F'\ttime: {end - start}s')

    # -----------------
    # resove duplicates
    # -----------------
    print('\n# -----------------\n# resolve duplicates\n# -----------------')
    start = time.time()

    model = resolve_duplicates(model,
                               check_reac = check_dupl_reac,
                               check_meta = check_dupl_meta,
                               remove_unused_meta = remove_unused_meta,
                               remove_dupl_reac = remove_dupl_reac,
                               replace_dupl_meta = remove_dupl_meta)

    end = time.time()
    print(F'\ttime: {end - start}s')

    # ----------
    # gapfilling
    # ----------

    print('\n# ----------\n# gapfilling\n# ----------')
    start = time.time()

    # construct a list of media
    # -------------------------

    media_list = []

    # load media from config file
    if media_path:
        media_list = read_media_config(media_path)

    # perform gapfilling
    # -----------------
        
    # ..........................................
    # @TODO 
    #   add an option for refinegems gapfilling
    #   separate option for cobra gapfilling 
    if len(media_list) > 0:
        # load universal model
        universal_model = load_model(universal,'cobra')
        # run gapfilling
        model = multiple_cobra_gapfill(model,universal_model,media_list,namespace,iterations=iterations, chunk_size=chunk_size, growth_threshold=growth_threshold)
    # ..........................................

    end = time.time()
    print(F'\ttime: {end - start}s')

    # ---------------------
    # dead ends and orphans
    # ---------------------
    # @TODO
    # currently no removal of dead ends and orphans as they may be interesting
    # for manual curation

    # ----------
    # save model
    # ----------
    name = F'{model.id}_clean'
    cobra.io.write_sbml_model(model, Path(dir,'step2-clean-up',name+'.xml'))

    # ---------------------------------
    # assess model quality using memote
    # ---------------------------------
        
    if memote:
        memote_path = Path(dir,'step2-clean-up',name+'.html')
        run_memote(model, 'html', return_res=False, save_res=memote_path, verbose=True)
