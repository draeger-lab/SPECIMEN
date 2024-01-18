"""Utility functions & Co. for handling cobra models.
"""

__author__ = 'Carolin Brune'

################################################################################
# requirements
################################################################################

import logging

from ..classes import reports

logger = logging.getLogger(__name__)

################################################################################
# variables
################################################################################

OBJECTIVE_FUNC_NAMES_GROWTH = ['Growth','GROWTH', 'growth',
                               'Biomass','BIOMASS', 'biomass']
MIN_GROWTH_RATE = 0.0005


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

# very similar to refinegems.biomass.test_biomass_presence(model)
# -> however that returns a list and not a string
# @DELETE, probably
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
