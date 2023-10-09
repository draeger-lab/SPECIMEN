"""Utility functions & Co. for handling cobra models.
"""

__author__ = 'Carolin Brune'

################################################################################
# requirements
################################################################################

import cobra
import logging
import os.path
from tqdm import tqdm
import warnings

from . import io
from ..classes import reports,medium

logger = logging.getLogger(__name__)

################################################################################
# variables
################################################################################

OBJECTIVE_FUNC_NAMES_GROWTH = ['Growth','GROWTH', 'growth',
                               'Biomass','BIOMASS', 'biomass']
MIN_GROWTH_RATE = 0.0005

# compartments
# ------------
# ....................
# @TODO
#     extension needed
# ....................
VALID_COMPARTMENTS = {'c': 'cytosol', 'e': 'extracellular space', 'p':'periplasm','y':'unknown compartment'}
COMP_MAPPING = {'c': 'c', 'e': 'e', 'p': 'p',
                'C_c': 'c', 'C_e': 'e', 'C_p': 'p',
                '':'y'}

# SBO terms
# ---------
sbo_bioch_terms = ["SBO:0000377", "SBO:0000399", "SBO:0000402", "SBO:0000403",
                   "SBO:0000660", "SBO:0000178", "SBO:0000200", "SBO:0000214",
                   "SBO:0000215", "SBO:0000217", "SBO:0000218", "SBO:0000219",
                   "SBO:0000220", "SBO:0000222", "SBO:0000223", "SBO:0000233",
                   "SBO:0000376", "SBO:0000401"]

################################################################################
# functions
################################################################################

def find_growth_obj_func(model):
    """Get the ID of the growth function of the model.

    Uses a list of standart names to search for a match and returns the first one.
    If no match is found, uses the objective expression instead.

    :param model: The model.
    :type model: cobra.Model

    :returns: The ID of the reaction identified/set as the growth.
    :rtype: string
    """
    # easy way: via name
    # get reaction names
    reacs = [_.id for _ in model.reactions]
    for name in OBJECTIVE_FUNC_NAMES_GROWTH:
        if name in reacs:
            logging.info(F'Growth function identified: {name}')
            return name
    # if nothing found, try via objective function
    logging.warning('No growth function found. Using the current objective function instead, but results might not be accurate.')
    return str(model.objective.expression).split(' ')[0].split('*')[1]


# @TODO
def reannotate_sbo_memote(model):
    """Remap the very specific SBO terms from SBOannotator to the high-level ones
    used by memote. Directly changes the model annotations.

    :param model: The model to be reannotated.
    :type model: cobra.Model

    :returns: The reannotated model
    :rtype: cobra.Model
    """

    # biochem reactions
    for r in model.reactions:
        if 'sbo' in r.annotation:
            if r.annotation['sbo'] in sbo_bioch_terms:
                r.annotation['sbo'] = "SBO:0000176"

    # @TODO: add transport reactions?

    return model


# dealing with compartments
# -------------------------

def are_compartment_names_valid(model):
    """Checks if the compartments of a model are valid.

    Note:
    Valid means their keys are part of VALID_COMPARTMENTS.

    :param model: A model.
    :type model: cobra.Model

    :returns: Result of the test. True if names are valid.
    :rtype: bool
    """

    for c in model.compartments.keys():
        if c not in VALID_COMPARTMENTS.keys():
            return False

    return True


def resolve_compartment_names(model):
    """Try to resolve invalid compartment names for a given model.

    Note:
    Valid means their keys are part of VALID_COMPARTMENTS.
    Valid mappings for solving the problem are in COMP_MAPPING.

    :param model: A model.
    :type model: cobra.Model

    :raises:  :class:`KeyError`: Unknown compartment detected. Cannot resolve problem.

    :returns: The updated model if possible
    :rtype: cobra.Model
    """

    # check if compartment names are valid
    if not are_compartment_names_valid(model):

        # check if mapping is possible
        if set(model.compartments.keys()).issubset(set(COMP_MAPPING.keys())):
            # for each metabolite rename the compartment
            for metabolite in model.metabolites:
                metabolite.compartment = COMP_MAPPING[metabolite.compartment]
            # add whole descriptions of the compartments to the model
            # note:
            #    only compartments IN the model will be added
            model.compartments = VALID_COMPARTMENTS

        else:
            raise KeyError(F'Unknown compartment {[_ for _ in model.compartments if _ not in COMP_MAPPING.keys()]} detected. Cannot resolve problem.')

    return model


# core-pan modelling
# ------------------
# ..............................................................................
# @TODO:
#     this part is still a working progress
#     currently heavily dependant on BiGG IDs as identifiers for reactions
# ..............................................................................

def get_intersecting_reactions(reacA, reacB, based_on='id'):
    """From two lists of reactions, get the intersection.
    If reactions are determined the same, the reaction object from the first
    reaction list is taken for the output list.

    The intersection can be done based on:
        'id' : use the reaction id for comparison

    :param reacA: List of reactions.
    :type reacA: list, list elements should be cobra.Reaction
    :param reacB: Another list of reactions.
    :type reacB: list, list elements should be cobra.Reaction
    :param based_on: Option on how to define reactions as equal.
        Default is 'id', all options can be found in the function description.
    :type based_on: string

    :raises:  :class:`ValueError`: Unknown option for intersection based_on.

    :returns: The intersection of the two lists.
    :rtype: list, list elements should be cobra.Reaction
    """

    match based_on:
        # match based on identifier
        case 'id':
            a = [_.id for _ in reacA]
            b = [_.id for _ in reacB]
            a_and_b = list(set(a) & set(b))
            return [_ for _ in reacA if _.id in a_and_b]

        case _:
            raise ValueError('Unknown option for intersection based_on: ',based_on)


def get_unique_reactions(reacA, reacB, based_on='id'):
    """From two lists of reactions, get the unique reactions of the first list.

    Unique reactions are determined by comparing:
        'id' : use the reaction id for comparison

    :param reacA: List of reactions.
    :type reacA: list, list elements should be cobra.Reaction
    :param reacB: Another list of reactions.
    :type reacB: list, list elements should be cobra.Reaction
    :param based_on: Option on how to define reactions as equal.
        Default is 'id', all options can be found in the function description.
    :type based_on: string

    :raises:  :class:`ValueError`: Unknown option for based_on.

    :returns: The reactions of list A that are not in list B.
    :rtype: list, list elements should be cobra.Reaction
    """

    match based_on:
        # match based on identifier
        case 'id':
            b = [_.id for _ in reacB]
            return [_ for _ in reacA if _.id not in b]

        case _:
            raise ValueError('Unknown option for based_on: ',based_on)


def remove_duplicate_reactions(reac, based_on='id'):
    """From a list of reactions, remove the duplicate reactions.

    Duplicate reactions are determined by comparing:
        'id' : use the reaction id for comparison

    :param reacA: List of reactions.
    :type reacA: list, list elements should be cobra.Reaction
    :param reacB: Another list of reactions.
    :type reacB: list, list elements should be cobra.Reaction
    :param based_on: Option on how to define reactions as equal.
        Default is 'id', all options can be found in the function description.
    :type based_on: string

    :raises:  :class:`ValueError`: Unknown option for based_on.

    :returns: The reactions of list A that are not in list B.
    :rtype: list, list elements should be cobra.Reaction
    """

    match based_on:
        # match based on identifier
        case 'id':
            reac_ids = list(set([_.id for _ in reac]))
            ids_entered = []
            reac_removed = []
            for r in reac:
                if r.id not in ids_entered:
                    reac_removed.append(r)
                    ids_entered.append(r.id)

            return reac_removed

        case _:
            raise ValueError('Unknown option for based_on: ',based_on)


def define_core_pan_reactions(reac_list, core, pan, based_on='id'):
    """Given three lists of new reactions, core reaction and pan reactions,
    sort the reactions of the first list into the latter two.

    Matching of the reactions can be based on:
        'id' : use the reaction id for comparison

    :param reac_list: The list of new reactions.
    :type reac_list: list, list elements should be cobra.Reaction
    :param core: The list of core reactions.
    :type core: list, list elements should be cobra.Reaction
    :param pan: The list of pan reactions.
    :type pan: list, list elements should be cobra.Reaction
    :param based_on: Option on how to define reactions as equal.
        Default is 'id', all options can be found in the function description.
    :type based_on: string

    :returns: The updated lists of core and pan reactions.
    :rtype: tuple, (core,pan)
    """


    if len(core) == 0:
        # nothing to compare => all is core
        core = reac_list
    else:
        temp = core
        # set core as reactions that can be found in the reac and the core list
        core = get_intersecting_reactions(reac_list,core, based_on)
        # extend pan by reactions no longer in core
        pan.extend(get_unique_reactions(temp,core,based_on))
        # extend pan by reactions from reac that are not in core
        pan.extend(get_unique_reactions(reac_list,core,based_on))
        # remove duplicates from pan
        pan = remove_duplicate_reactions(pan,based_on)

    return (core,pan)


def build_core_pan_model(file_list, based_on='id', remove_genes=True):
    """From a list of file paths to models, create a core-pan model from their reactions.

    The core-pan model is created by merging all reactions of the models into
    one model. Reactions that are found in all models are characterised as 'core',
    while the others are characterised as 'pan'. This characterisation is saved
    in the notes-field unter 'core-pan'.

    Note:
        To make the final model 'cleaner', each model is checked for their
        compartment names using `resolve_compartment_names()`.

    Matching of the reactions can be based on:
        'id' : use the reaction id for comparison

    :param file_list: List of file paths to GEMs.
    :type file_list: list
    :param based_on: Option on how to define reactions as equal.
        Default is 'id', all options can be found in the function description.
    :type based_on: string
    :param remove_genes: Option to remove the genes / gene reaction rules from the model.
        Default is True, as the core-pan differentiation is based on ID and therefore on the reactions.
    :type remove_genes: bool

    :returns: The core-pan model.
    :rtype: cobra.Model
    """

    # step 1: define the core and pan reactions
    #         for a given set of models
    core = []
    pan = []
    for f in tqdm(file_list):
        # load model
        model = io.read_model_cobra(f)
        # resolve compartment naming issues
        model = resolve_compartment_names(model)
        # idenfy core-pan reactions
        reac = model.reactions
        core,pan = define_core_pan_reactions(reac,core,pan, based_on)

    # step 2: label the reactions accordingly
    # label with core / pan
    for reac in core:
        reac.notes['pan-core'] = 'core'
    for reac in pan:
        reac.notes['pan-core'] = 'pan'

    # step 3: build a core-pan model from the reactions
    cp_model = cobra.Model('core_pan_model')
    cp_model.add_reactions(core)
    cp_model.add_reactions(pan)

    # step 4: remove genes (optional)
    if remove_genes:
        cobra.manipulation.delete.remove_genes(cp_model, cp_model.genes, remove_reactions=False)

    return cp_model


# core_pan_analysis
# -----------------

def compare_to_core_pan(model, cp_model, based_on='id'):
    """Perform a pan-core model analysis by comparing an input model and a pan-core model.

    Currently only checks the reactions. Depends highly on the pan-core model being
    build with `build_core_pan_model()`.

    Matching/Comparison can be based on:
        'id' : use the reaction id for comparison

    @TODO
        - implements more options for comparison
        - include checking for metabolites

    :param model: The model to perform the analysis on.
    :type model: cobra.Model
    :param cp_model: The pan-core model to compare to.
    :type cp_model: cobra.Model
    :param based_on: Option on how to define reactions as equal.
        Default is 'id', all options can be found in the function description.
    :type based_on: string

    :returns: A report of the comparison.
    :rtype: report.PanCoreAnalysisReport
    """

    results = reports.PanCoreAnalysisReport(model)

    # reactions
    # ---------

    # separate cp_model reactions into core and pan list
    core_reac_list = [_ for _ in cp_model.reactions if _.notes['pan-core']=='core']
    pan_reac_list = [_ for _ in cp_model.reactions if _.notes['pan-core']=='pan']

    # compare model to the core and pan reaction list
    results.core_reac = get_intersecting_reactions(model.reactions, core_reac_list, based_on)
    results.pan_reac = get_intersecting_reactions(model.reactions, pan_reac_list, based_on)
    results.novel_reac = get_unique_reactions(model.reactions, cp_model.reactions, based_on)

    # ...................
    # metabolites
    # -----------
    # @TODO
    # ...................

    return results


# statistical analysis
# --------------------
def generate_statistics(model):

    # reactions
    # ---------

    # total number of reactions
    total_reac = len(model.reactions)

    # number of reactions from template, KEGG / MetaNetX, gapfilling
    reac_origin_counts = {'via template':0, 'via MetaNetX':0, 'via KEGG':0, 'via gapfilling':0, 'else':0}
    reac_with_gpr = 0
    for reac in model.reactions:
        # get origin of reaction (based on workflow notation)
        if 'creation' in reac.notes.keys():
            if reac.notes['creation'] in reac_origin_counts.keys():
                reac_origin_counts[reac.notes['creation']] += 1
            else:
                reac_origin_counts['else'] += 1
        else:
            reac_origin_counts['else'] += 1

        # check for GPR
        if len(reac.genes) > 0:
            reac_with_gpr += 1

    # metabolites
    # -----------

    # total number of reactions
    total_meta = len(model.metabolites)

    # genes
    # -----

    # total number of genes
    total_gene = len(model.genes)

    return reports.ModelStatisticsReport(model.id, total_reac, total_meta, total_gene, reac_origin_counts, reac_with_gpr)
