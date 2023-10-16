"""Perform part 3 of the refinement: clean-up.
"""

__author__ = 'Carolin Brune'

################################################################################
# requirements
################################################################################

import time
from pathlib import Path
import warnings
import subprocess

import pandas as pd
import numpy as np
import math
import cobra
import copy
import re
import sys

from ...classes import medium
from ... import util
# further required programs:
#        - MEMOTE,  tested with version 0.13.0

################################################################################
# variables
################################################################################

NH_PATTERN = re.compile('nh[3-4]')

################################################################################
# functions
################################################################################

def check_direction(model,data_file):
    """Check the direction of newly created reactions (01-extention) by searching for matching MetaCyc,
    KEGG and MetaNetX IDs as well as EC number in a downloaded BioCyc (MetaCyc)
    database table (need to contain at least the following columns:
    Reactions (MetaCyc ID),EC-Number,KEGG reaction,METANETX,Reaction-Direction.

    note: checks only creations that do not contain the notes['creation'] == 'via template'

    :param model:      The GEM containing the reactions to be check for direction.
    :type  model:      cobra.Model
    :param data_file:  Path to the MetabCyc (BioCyc) smart table.
    :type  data_file:  string
    :returns:          The updated model.
    :rtype:            cobra.Model
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


def complete_BioMetaCyc(model):
    """Check for existing MetaCyc and BioCyc annotations for metabolites and
    reactions and generate them from the other if one of the two is missing.

    :param model: A genome-scale model to be checked for complete BioCyc/MetaCyc annotations.
    :type  model: cobra.Model
    :returns:     The updated model.
    :rtype:       cobra.Model
    """

    # reactions
    # ---------
    for reac in model.reactions:
        # case 1: MetaCyc but not BioCyc
        if 'metacyc.reaction' in reac.annotation and not 'biocyc' in reac.annotation:
            # add database (META) information to get BioCyc ID
            if isinstance(reac.annotation['metacyc.reaction'], list):
                reac.annotation['biocyc'] = ['META:' + _ for _ in reac.annotation['metacyc.reaction']]
            else:
                reac.annotation['biocyc'] = 'META:' + reac.annotation['metacyc.reaction']
        # case 2: BioCyc, but no MetaCyc
        elif 'biocyc' in reac.annotation and not 'metacyc.reaction' in reac.annotation:
            # if there are multiple metacyc.reaction annotation
            if isinstance(reac.annotation['biocyc'], list):
                add_anno = []
                for biocyc_anno in reac.annotation['biocyc']:
                    if ':' in biocyc_anno:
                        # exclude organism identifier from ID to get MetaCyc ID
                        add_anno.append(biocyc_anno.split(':')[1])
                    else:
                        # high possibility that information is faulty - do not use it
                        print(F'\n\nWarning: Unexpected BioCyc annotation {biocyc_anno} for reaction {reac.id}')
                reac.annotation['metacyc.reaction'] = add_anno
            # if there is only one
            else:
                if ':' in reac.annotation['biocyc']:
                    # exclude organism identifier from ID to get MetaCyc ID
                    reac.annotation['metacyc.reaction'] = reac.annotation['biocyc'].split(':')[1]
                else:
                    # high possibility that information is faulty - do not use it
                    print(F'\n\nWarning: Unexpected BioCyc annotation {reac.annotation["biocyc"]} for reaction {reac.id}')
        # case 3: both or no = skip
        else:
            continue

    # metabolites
    # -----------
    for meta in model.metabolites:
        # case 1: MetaCyc but not BioCyc
        if 'metacyc.compound' in meta.annotation and not 'biocyc' in meta.annotation:
            # add database (META) information to get BioCyc ID
            if isinstance(meta.annotation['metacyc.compound'],list):
                meta.annotation['biocyc'] = ['META:' + _ for _ in meta.annotation['metacyc.compound']]
            else:
                meta.annotation['biocyc'] = 'META:' + meta.annotation['metacyc.compound']
        # case 2: BioCyc, but no MetaCyc
        elif 'biocyc' in meta.annotation and not 'metacyc.compound' in meta.annotation:
            # if there are multiple metacyc.compound annotations
            if isinstance(meta.annotation['biocyc'], list):
                add_anno = []
                for biocyc_anno in meta.annotation['biocyc']:
                    if ':' in biocyc_anno:
                        # exclude organism identifier from ID to get MetaCyc ID
                        add_anno.append(biocyc_anno.split(':')[1])
                    else:
                        # high possibility that information is faulty - do not use it
                        print(F'\n\nWarning: Unexpected BioCyc annotation {biocyc_anno} for metabolite {meta.id}')
                meta.annotation['metacyc.compound'] = add_anno
            # if there is only one
            else:
                if ':' in meta.annotation['biocyc']:
                    # exclude organism identifier from ID to get MetaCyc ID
                    meta.annotation['metacyc.compound'] = meta.annotation['biocyc'].split(':')[1]
                else:
                    # high possibility that information is faulty - do not use it
                    print(F'\n\nWarning: Unexpected BioCyc annotation {meta.annotation["biocyc"]} for metabolite {meta.id}')
        # case 3: both or no = skip
        else:
            continue

    return model


# handling duplicates
# -------------------
def resolve_duplicate_reactions(model, based_on='reaction', remove_reac = True):
    """Resolve and remove duplicate reaction based on their reaction equation
    and matching database identifiers. Only if all match or are nan will one of
    the reactions be removed.

    :param model:       The model to be checked for duplicate reactions.
    :type  model:       cobra.Model
    :param based_on:    Starting point of the comparison for duplicates. Default is 'reaction'. Currently not advisable to change this parameter.
    :type  based_on:    string
    :param remove_reac: Option to keep or remove found duplicates.
    :type  remove_reac: bool, default True.
    :returns:           The updated model.
    :rtype:             cobra.Model
    """

    # get annotation and compartment information
    anno_reac = []
    for r in model.reactions:
        anno_reac.append({'id':r.id, 'compartment':str(r.compartments), 'reaction':r.reaction} | r.annotation)
    df_reac = pd.DataFrame.from_dict(anno_reac)

    # check if based_on is valid
    if not based_on in df_reac.columns.tolist():
        warnings.warn(F'Warning: Annotation column {based_on} does not exists. Search for duplicates will be skipped.')
        return model

    # set basic parameters
    skip_cols = ['id','compartment','bigg.reaction','reaction',based_on]
    colnames = df_reac.columns.tolist()

    for c in df_reac.groupby('compartment'):
        # note: using groupby drops nans
        for mnx in c[1].groupby(based_on):
            # find possible duplicates
            dupl = True
            annotations = {}
            if len(mnx[1]) > 1:
                # check annotations
                for col in [_ for _ in colnames if not _ in skip_cols]:
                    if len(mnx[1][col].dropna().value_counts()) < 2:
                        annotations[col] = mnx[1][col].dropna().explode().unique().tolist()
                    else:
                        dupl = False
                        break

                # if duplicate found
                if dupl:
                    # choose reaction to keep
                    keep_reac = model.reactions.get_by_id(mnx[1]['id'].tolist()[0])

                    # resolve annotations
                    for key, value in annotations.items():
                        if len(value) > 0 and not key in keep_reac.annotation:
                            keep_reac.annotation[key] = value

                    # combine gene reaction rules
                    for r_id in mnx[1]['id'].tolist()[1:]:
                        keep_reac.gene_reaction_rule = keep_reac.gene_reaction_rule + model.reactions.get_by_id(r_id).gene_reaction_rule
                        if remove_reac:
                            model.reactions.get_by_id(r_id).delete()
                            print(F'\tDuplicate reaction {r_id} found. Combined to {keep_reac.id} and deleted.')
                        else:
                            print(F'\tDuplicate reaction {r_id} found. Combined to {keep_reac.id}. Deletion skipped.')
    return model


def resolve_duplicate_metabolites(model, based_on='metanetx.chemical', replace=True):
    """Resolve duplicate metabolites in a model. Metabolites are considered
    duplicate if they share the same annotations (same or nan).
    Note: Depending on the starting database, the results might differ.

    :param model:    The model to be checked for duplicate metabolites.
    :type  model:    cobra.Model
    :param based_on: Starting point of the comparison.
    :param replace:  Option to replace the duplicate metabolites or only report them.
    :type  replace:  bool, default is True
    :returns:        The updated model.
    :rtype:          cobra.Model
    """

    # get annotation and compartment information
    anno_meta = []
    for m in model.metabolites:
        anno_meta.append({'id':m.id, 'compartment':m.compartment} | m.annotation)
    df_meta = pd.DataFrame.from_dict(anno_meta)

    # check if based_on is valid
    if not based_on in df_meta.columns.tolist():
        warnings.warn(F'Warning: Annotation {based_on} not found. Search for metabolite duplicates skipped.')
        return model

    # set basic parameters
    skip_cols = ['id','compartment','bigg.metabolite',based_on]
    colnames = df_meta.columns.tolist()
    objective_function = util.cobra_models.find_growth_obj_func(model)

    for c in df_meta.groupby('compartment'):
        # note: using groupby drops nans
        # note: mnx as starting point was choosen, as it seems to have the most annotations (easy to get)
        for mnx in c[1].groupby(based_on):

            # find possible duplicates
            dupl = True
            annotations = {}
            if len(mnx[1]) > 1:
                for col in [_ for _ in colnames if not _ in skip_cols]:
                    # check annotation, if current group truly consists of duplicates
                    if len(mnx[1][col].dropna().value_counts()) < 2:
                        annotations[col] = mnx[1][col].dropna().explode().unique().tolist()
                    else:
                        dupl = False
                        break
            # no duplicates = check next
            else:
                continue

            # if duplicate is found:
            if dupl:
                # either replace ...
                if replace:
                    # choose metabolite to keep
                    keep_meta = model.metabolites.get_by_id(mnx[1]['id'].tolist()[0])
                    # resolve annotations
                    for key, value in annotations.items():
                        if len(value) > 0 and not key in keep_meta.annotation:
                            keep_meta.annotations[key] = value
                    # note: charge and formula should be in valid range and be corrected by MCC (if needed)

                    # replace duplicates with the metabolite to be kept
                    #     to ensure consistency, only delete duplicate metabolites, which
                    #     do NOT share ANY reactions

                    # retrieve reactions for metabolites set for keeping
                    keep_reac = [_.id for _ in model.metabolites.get_by_id(keep_meta.id).reactions]
                    # iterate over metabolites set for deletion
                    for del_meta_id in mnx[1]['id'].tolist()[1:]:
                        # retrieve reaction for metabolites set for deletion
                        del_reac = [_.id for _ in model.metabolites.get_by_id(del_meta_id).reactions]
                        # get intersection of reactions (keep + del)
                        reac_intersec = list(set(keep_reac) & set(del_reac))

                        # if intersection empty, metabolite is with a high probability indeed a duplicate
                        # Special case: NH3 / NH4
                        #    intersection does not have to be emtpy, know 'problem' caused by CarveMe
                        if len(reac_intersec) == 0 or all([re.search(NH_PATTERN,_) for _ in [keep_meta.id,del_meta_id]]):
                            # automated deletion is only advisable, if consistency can be retained
                            perform_deletion = True
                            with model as model_del:

                                # if the special case is detected ...
                                if all([re.search(NH_PATTERN,_) for _ in [keep_meta.id,del_meta_id]]):
                                    print(F'\tSpecial case -Duplicate NH4/NH3- detected.\n\tTrying to solve by additionally removing reactions containing both metabolites.')
                                    # ... remove reactions with nh3 and nh4 both present
                                    for del_reac_id in reac_intersec:
                                        # if objective_function is part of the set
                                        # automated deletion is (currently) not possible
                                        if del_reac_id == objective_function:
                                            perform_deletion = False
                                            break
                                        model_del.reactions.get_by_id(del_reac_id).remove_from_model()
                                    # set the metabolites to be deleted to be the one NOT in the objective functions
                                    # to avoid inconsistencies
                                    if del_meta_id in [_.id for _ in model_del.reactions.get_by_id(objective_function).metabolites]:
                                        temp = del_meta_id
                                        del_meta_id = keep_meta.id
                                        keep_meta = model_del.metabolites.get_by_id(temp)

                                # try replacing metabolite with the kept duplicate ...
                                for reac in model_del.metabolites.get_by_id(del_meta_id).reactions:

                                    reac.reaction = reac.reaction.replace(del_meta_id,keep_meta.id)

                                    # skip if objective function is found
                                    if reac.id == objective_function:
                                        continue

                                    # check if consistency is still intact
                                    balance_test = reac.check_mass_balance()
                                    if not reac.boundary and len(balance_test) > 0:
                                        # try fixing H-balance
                                        if 'H' in balance_test.keys():
                                            # ..............................
                                            # TODO:
                                            #    get H according to compartment
                                            #    current implementation relies heavily
                                            #    on 'correct' use input: compartment should have format C_c or c (C_p, p, C_e, e etc.)
                                            # ..............................
                                            reac_comp = reac.compartments.pop()[-1]
                                            if reac_comp == 'c':
                                                reac.subtract_metabolites({'h_c':balance_test['H']})
                                            elif reac_comp == 'p':
                                                reac.subtract_metabolites({'h_p':balance_test['H']})
                                            elif reac_comp == 'e':
                                                reac.subtract_metabolites({'h_e':balance_test['H']})
                                            else:
                                                perform_deletion = False
                                                break
                                        # ..............................
                                        # TODO:
                                        #    fix other possible problems
                                        # ..............................

                                        # finally, check balance again (continue only if fixed, else break)
                                        if len(reac.check_mass_balance()) > 0:
                                            perform_deletion = False
                                            break

                                    else:
                                        continue

                                # if not problems are found, duplicate is removed
                                if perform_deletion:
                                    model = model_del.copy()
                                    print(F'\tDuplicate metabolite {del_meta_id} found. Replaced with {keep_meta.id}.')
                                # if problems are not solvable, duplicate is kept and only reported
                                else:
                                    print(F'\tDuplicate metabolite {del_meta_id} found (duplicate to {keep_meta.id} based on annotation).\n\t\tAutomated deletion not possible due to problems with consistency.')

                        # else, metabolite is kept
                        #       since it might be an isomer, elongation, or other explanation
                        #       for the same annotation
                        else:
                            print(F'\tDuplicate metabolite {del_meta_id} found (duplicate to {keep_meta.id} based on annotation).\n\t\tKept, as reaction containing both metabolites was found.')


                # ... or only report duplicates
                else:
                    print(F'\tDuplicate metabolite(s) {", ".join(mnx[1]["id"].tolist())} found.')


    return model


def resolve_duplicates(model, check_reac = True, check_meta = 'default', replace_dupl_meta=True, remove_unused_meta = False, remove_dupl_reac = True):
    """Resolve and remove (optional) duplicate metabolites and reactions in the model.

    :param model:              The model to be checked for duplicates.
    :type  model:              cobra.Model
    :param check_reac:         Option to check for duplicate reactions.
    :type  check_reac:         bool, default is True (check reactions).
    :param check_meta:         Option to check for duplicate metabolites.
    :type  check_meta:         bool, default is True (check metabolites).
    :param replace_dupl_meta:  Option to replace duplicate metabolites.
    :type  replace_dupl_meta:  bool, default is True (Replace metabolites with duplicate).
    :param remove_unused_meta: Option to remove unused metabolites (metabolites that are not taking part in any reaction).
    :type  remove_unused_meta: bool, default is False (remove unused metabolites).
    :param remove_dupl_reac:   Option to remove duplicate reactions.
    :type  remove_dupl_reac:   bool, default is True (remove duplicate reactions).
    :returns:                  The updated model.
    :rtype:                    cobra.Model
    """

    # resolve duplicate metabolites
    if check_meta == 'default':
        # resolve duplicates starting with the metanetx.chemical database identifiers
        model = resolve_duplicate_metabolites(model, replace=replace_dupl_meta)
    elif check_meta == 'exhaustive':
        # resolve duplicates by starting at every database identifer one after another
        # note: bigg and sbo are skipped as sbo gives not much information and bigg is
        #       usually the one that differs (naming issue)
        for colname in [_ for _ in df_meta.columns.tolist() if not _ in ['id','compartment','bigg.metabolite','sbo']]:
            model = resolve_duplicate_metabolites(model,colname,replace=replace_dupl_meta)
    elif check_meta == 'skip':
        print('\tSkip check for duplicate metabolites.')
    else:
         warnings.warn(F'Warning: Unknown option for metabolites duplicate checking {check_meta}. Search for metabolite duplicates skipped.')

    # remove now unused metabolites
    if remove_unused_meta:
        model,removed = cobra.manipulation.delete.prune_unused_metabolites(model)
        print(F'\tThe following metabolites () have been removed: {", ".join([x.id for x in removed])}')

    # resolve duplicate reactions
    if check_reac:
        model = resolve_duplicate_reactions(model, based_on='reaction', remove_reac =remove_dupl_reac)

    return model


# gap filling
# -----------

def single_cobra_gapfill(model, universal, medium, growth_threshold = 0.05):
    """Attempt gapfilling (with cobra) for a given model to allow growth on a given
    medium.

    :param model: The model to perform gapfilling on.
    :type model: cobra.Model
    :param universal: A model with reactions to be potentially used for the gapfilling.
    :type universal: cobra.model
    :param medium: A medium the model should grow on.
    :type medium: medium.Medium
    :param growth_threshold: Minimal rate for the model to be considered growing.
    :type: float

    :returns: List of reactions to be added to the model to allow growth
        or True, if the model already grows.
    :rtype: list or bool
    """
    # perform the gapfilling
    solution = []
    with model as model_copy:
        # set medium model should grow on
        medium.validate_for_model(model)
        model_copy.medium = medium.export_to_cobra()
        # if model does not show growth (depending on threshold), perform gapfilling
        if model_copy.slim_optimize() < growth_threshold:
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


def cobra_gapfill_wrapper(model, universal, medium, iterations=3, chunk_size=10000, growth_threshold = 0.05):
    """Wrapper for single_cobra_gapfill().

    Either use the full set of reactions in universal model by setting iteration to
    0 or None or use them in randomized chunks for faster runtime (especially useful
    on laptops). Note: when using the second option, be aware that this does not
    test all reaction combinations exhaustively (heuristic approach!!!).

    :param model: The model to perform gapfilling on.
    :type model: cobra.Model
    :param universal: A model with reactions to be potentially used for the gapfilling.
    :type universal: cobra.model
    :param medium: A medium the model should grow on.
    :type medium: medium.Medium
    :param iterations: Number of iterations for the heuristic version.
        Default is 3, is 0 or None is given, uses full set of reactions.
    :type iterations: int
    :param chunk_size: Number of reactions to be used for gapfilling at the same time.
        Default is 10000. If None or 0 is given, use full set, not heuristic.
    :type chunk_size: int
    :param growth_threshold: Minimal rate for the model to be considered growing.
    :type: float

    :returns: The updated (gapfilled) model.
    :rytpe: cobra.Model
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
                solution = single_cobra_gapfill(model, subset_reac, medium, growth_threshold)

                if (isinstance(solution,bool) and solution) or (isinstance(solution, list) and len(solution) > 0):
                    break

            if (isinstance(solution,bool) and solution) or (isinstance(solution, list) and len(solution) > 0):
                break

    # use the whole reactions content of the universal model at once
    #     not advised for Laptops and PCs with small computing power
    #     may take a long time
    else:
        solution = single_cobra_gapfill(model, universal, medium, growth_threshold)

    # if solution is found add new reactions to model
    if isinstance(solution, list) and len(solution) > 0:
        for reac in solution[0]:
            reac.notes['creation'] = 'via gapfilling'
        print(F'Adding {len(solution[0])} reactions to model to ensure growth on medium {medium.name}.')
        model.add_reactions(solution[0])

    return model


def multiple_cobra_gapfill(model, universal, media_list, growth_threshold = 0.05, iterations=3, chunk_size=10000):
    """Perform single_cobra_gapfill() on a list of models.

    :param model: The model to perform gapfilling on.
    :type model: cobra.Model
    :param universal: A model with reactions to be potentially used for the gapfilling.
    :type universal: cobra.model
    :param media: A list of media the model should grow on.
    :type media: list of medium.Medium
    :param growth_threshold: Minimal rate for the model to be considered growing.
    :type: float
    :param iterations: Number of iterations for the heuristic version.
        Default is 3, is 0 or None is given, uses full set of reactions.
    :type iterations: int
    :param chunk_size: Number of reactions to be used for gapfilling at the same time.
        Default is 10000. If None or 0 is given, use full set, not heuristic.
    :type chunk_size: int

    :returns: The gapfilled model, of possible.
    :rtype: cobra.model
    """

    for medium in media_list:
        model = cobra_gapfill_wrapper(model,universal,medium, iterations, chunk_size, growth_threshold)

    return model


# run
# ---

def run(model, dir, biocyc_db=None, check_dupl_reac = False,
        check_dupl_meta = 'default',
        remove_unused_meta = False, remove_dupl_reac = False, remove_dupl_meta = False,
        universal = None, media_db = None, load_media = None, external_media = None,
        change_to_aerobic = None, change_to_anaerobic = None,
        add_casamino = None, growth_threshold = 0.05,
        iterations=3, chunk_size=10000,
        memote = False):

    """Perform the second refinement step, cleanup, on a model.

    The second refinement step resolves the following issues:
    - (optional) checking direction of reactions with BioCyc
    - complete BioCyc/MetaCyc annotation inconsistencies
    - find and/or resolve duplicates (reactions and metabolites)
    - gapfilling using cobra

    :param model: The Path to an sbml model.
    :type model: string
    :param dir: Path to the directory of the output.
    :type dir: string
    :param biocyc_db: Path to the BioCyc/MetaCyc reaction database file.
        Default is None, which leads to skipping the direction check.
    :type biocyc_db: string, optional
    :param check_dupl_reac: Option to check for duplicate reactions.
        Default is False
    :type check_dupl_reac: bool
    :param check_dupl_meta: Option to check for duplicate metabolites.
        Default is 'default' (check based on MetaNetX first).
        Further options include 'skip' and 'exhaustive' (check for all possibilities).
    :type check_dupl_meta: bool
    :param remove_unused_meta: Option to remove unused metabolites from the model.
        Default is False.
    :type remove_unused_meta: bool
    :param remove_dupl_reac: Option to remove the duplicate reactions.
        Deafult is False. True is only applied, if check_dupl_reac is also True.
    :type remove_dupl_reac: bool
    :param remove_dupl_meta: Option to remove the duplicate metabolites.
        Deafult is False. True is only applied, if check_dupl_meta is also True.
    :type remove_dupl_meta: bool

    @TODO missing docstrings

    :param iterations: Number of iterations for the heuristic version.
        Default is 3, is 0 or None is given, uses full set of reactions.
    :type iterations: int
    :param chunk_size: Number of reactions to be used for gapfilling at the same time.
        Default is 10000. If None or 0 is given, use full set, not heuristic.
    :type chunk_size: int

    :param memote: Option to run memote on the cleaned model.
        Default is False. Requires Memote to be installed and runnable from the command line.
    :type memote: bool
    """

    # -------------------
    # checking parameters
    # -------------------
    if not check_dupl_meta in ['default','skip','exhaustive']:
        raise ValueError('Unknown option {check_dupl_meta} for checking duplicate metabolite. Use one of: default, skip, exhaustive')
        sys.exit(1)

    # -------------
    # start program
    # -------------
    print('\nrefinement step 2: clean-up\n################################################################################\n')

    # -----------------------
    # create output directory
    # -----------------------

    # make sure given directory path ends with '/'
    if not dir.endswith('/'):
        dir = dir + '/'

    try:
        Path(F"{dir}step2-clean-up/").mkdir(parents=True, exist_ok=False)
        print(F'Creating new directory {F"{dir}step2-clean-up/"}')
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

    # load media from database
    if load_media:
        db = medium.load_media_db(media_db)

        # load medium as it is
        if load_media:
            for medium_name in load_media:
                if medium_name in db.keys():
                    media_list.append(db[medium_name])
                else:
                    raise KeyError('Medium not found in database: ', medium_name)

        # load medium + make aerobic
        if change_to_aerobic:
            for medium_name in change_to_aerobic:
                new_medium = copy.deepcopy(db[medium_name])
                new_medium.make_aerobic()
                new_medium.name = new_medium.name + ' (aerobic)'
                media_list.append(new_medium)

        # load medium + make anaerobic
        if change_to_anaerobic:
            for medium_name in change_to_anaerobic:
                new_medium = copy.deepcopy(db[medium_name])
                new_medium.make_anaerobic()
                new_medium.name = new_medium.name + ' (anaerobic)'
                media_list.append(new_medium)

        # load medium + add CASAMINO_ACIDS
        if add_casamino:
            for medium_name in add_casamino:
                media_list.append(db[medium_name] + medium.CASAMINO_ACIDS)

    # load external media
    if external_media:
        for filepath in external_media:
            tmp_table = pd.read_csv(filepath, sep=';', header=0)
            media_list.append(medium.from_table(tmp_table))

    # perform gapfilling
    # -----------------
    if len(media_list) > 0:
        # load universal model
        universal_model = util.io.read_model_cobra(universal)
        # run gapfilling
        model = multiple_cobra_gapfill(model,universal_model,media_list,iterations=iterations, chunk_size=chunk_size, growth_threshold=growth_threshold)

    end = time.time()
    print(F'\ttime: {end - start}s')

    # ---------------------
    # dead ends and orphans
    # ---------------------
    # no removal of dead ends and orphans as they may be interesting
    # for manual curation

    # ----------
    # save model
    # ----------
    name = F'{model.id}_clean'
    cobra.io.write_sbml_model(model, F'{dir}step2-clean-up/{name}.xml')

    # ---------------------------------
    # assess model quality using memote
    # ---------------------------------

    if memote:
        print('\n# -------------------\n# analyse with MEMOTE\n# -------------------')
        start = time.time()
        draft_path = F'{dir}step2-clean-up/{name}.xml'.replace(" ", "\ ")
        memote_path = F'{dir}step2-clean-up/{name}.html'.replace(" ", "\ ")
        subprocess.run([F'memote report snapshot --filename {memote_path} {draft_path}'], shell=True)
        end = time.time()
        print(F'\ttotal time: {end - start}s')
