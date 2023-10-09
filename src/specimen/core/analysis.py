"""Analyse a model (part 5 of the workflow).
"""

__author__ = 'Carolin Brune'

################################################################################
# requirements
################################################################################

import cobra
import copy
import re
import time
import warnings

from pathlib import Path

from .. import util, classes

################################################################################
# functions
################################################################################

# growth_analysis
# ---------------
def growth_simulation(model, media_list):
    """Perform growth simulation for a list of media.

    :param model: The model for simulation.
    :type model: cobra.Model
    :param media_list: List of media to be tested.
    :type media_list: list of classes.medium.Medium

    :returns: The results.
    :rtype: dict, classes.medium.name:growth_rate
    """

    results = {}
    with model as m:
        for medium in media_list:
            medium.validate_for_model(m)
            m.medium = medium.export_to_cobra()
            results[medium.name] = m.slim_optimize()

    return results


def auxotrophies_simulation(model, media_list, minimal=None, show_names=False):
    """Perform auxotrophy test for amino acids on a model and a list of media.

    Tests, if the model growths on the media and if and with what fluxes the
    20 proteinogenic amino acids are produced by temporarily adding a
    sink reaction for each of the amino acids to the model as the objective function.

    :param model: The model to be tested.
    :type model: cobra.Model
    :param media_list: list of media to be tested.
    :type media_list: list, elements of type classes.medium.Medium
    :param minimal: A previously generated minimal medium.
    :type minimal: classes.medium.Medium, optional
    :param show_names: Use the names of the amino acids as axis labels instead of the BiGG IDs.
    :type show_names: bool

    :returns: The results of the tests.
    :rtype: nested dict, {<medium.name>:{<aa.name>:growth_rate}}
    """

    results = {}

    if minimal:
        media_list.append(minimal)

    # test the auxotrophies for all media
    old_medium = model.medium
    for med in media_list:
        auxotrophies = {}
        med.validate_for_model(model)
        # then iterate over all amino acids
        for aa in classes.medium.AMINO_ACIDS.compounds.keys():
            with model as m:
                # first set the medium
                m.medium = med.export_to_cobra()
                # check if amino acids is in model
                if not aa+'_c' in [_.id for _ in m.metabolites]:
                    warnings.warn(F'The following amino acid {aa} could not be found in the cytosol.\nIdentifier {aa+"c"} not found in model.metabolites.')
                    continue
                # else
                # create a pseudo reaction -> a sink reaction for the amino acid
                # to use as the new objective
                m.add_boundary(m.metabolites.get_by_id(F'{aa}_c'), type='sink', reaction_id=F'sink_{aa}_c_tmp')
                m.objective = F'sink_{aa}_c_tmp'
                # if existent, close the exchange reaction
                if F'EX_{aa}_e' in [_.id for _ in m.exchanges]:
                    m.reactions.get_by_id(F'EX_{aa}_e').lower_bound = 0.0
                    m.reactions.get_by_id(F'EX_{aa}_e').upper_bound = 0.0
                # and calculate the new objective
                if show_names:
                    auxotrophies[classes.medium.AMINO_ACIDS.compounds[aa].name] = m.slim_optimize()
                else:
                    auxotrophies[aa] = m.slim_optimize()

            # add the current test results to the list of all results
            results[med.name] = auxotrophies

    model.medium = old_medium

    return results


def analyse_growth(model, db_path=None, load_media=None, change_to_aerobic=None,
                   change_to_anaerobic=None, add_casamino=None, external_media=None,
                   test_min_medium=False, growth_rate=0.5,
                   test_aa_auxotrophies=True):
    """Analyse the growth of a model on different media.

    :param model: The model to be tested.
    :type model: cobra.Model
    :param db_path: Path to a media database file.
    :type db_path: string (optional)
    :param load_media: List of medium identifiers to be loaded from the database.
    :type load_media: list, elements are strings (optional)
    :param change_to_aerobic: List of medium identifiers to be loaded from the database
        and changed to fullfill the aerobic criteria.
    :type change_to_anaerobic: list, elements are strings (optional)
    :param change_to_aerobic: List of medium identifiers to be loaded from the database
        and changed to fullfill the anaerobic criteria.
    :type change_to_anaerobic: list, elements are strings (optional)
    :param add_casamino: List of medium identifiers to be loaded from the database
        Afterwards, add the casamino acids to the medium.
    :type add_casamino: list, elements are strings (optional)
    :param external_media: List of paths to files containing exactly one medium description.
        Offers the user the possibility to use personal media without adding them into a database.
    :type external_media: list, elements are strings (optional)
    :param test_min_medium: Option to search for a minimal medium based on the
        exchange reactions. This can take a while,if the number of exchange
        reactions is large. Default is False. Requires the :growth_rate:
        parameter to be higher than 0.0.
    :type test_min_medium: bool
    :param growth_rate: Growth rate that should be reached by the minimal medium.
        Default is 0.5. Rule of thumb: higher values lead to faster runtime.
    :type growth_rate: float, should be > 0.0.
    :param test_aa_auxotrophies: Enable the auxotrophy test.
        Default is True.
    :type test_aa_auxotrophies: bool

    :returns: The report of the analysis
    :rtype: classes.reports.GrowthAnalysisReport
    """


    if (not db_path and not load_media) and not external_media:
        warnings.warn('Found no media to test growth on, returning None')
        return None

    # construct a list of media
    # -------------------------

    media_list = []

    # load media from database
    if load_media:
        db = classes.medium.load_media_db(db_path)

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
                new_medium.name = new_medium.name + '_aerobic'
                media_list.append(new_medium)

        # load medium + make anaerobic
        if change_to_anaerobic:
            for medium_name in change_to_anaerobic:
                new_medium = copy.deepcopy(db[medium_name])
                new_medium.make_anaerobic()
                new_medium.name = new_medium.name + '_anaerobic'
                media_list.append(new_medium)

        # load medium + add CASAMINO_ACIDS
        if add_casamino:
            for medium_name in add_casamino:
                media_list.append(db[medium_name] + classes.medium.CASAMINO_ACIDS)

    # load external media
    if external_media:
        for filepath in external_media:
            tmp_table = pd.read_csv(path, sep=';', header=0)
            media_list.append(from_table(tmp_table))

    # performs tests
    # --------------
    # test the growth
    growth_sim_results = growth_simulation(model, media_list)

    # get minimal medium based on exchanges
    if test_min_medium:
        minimal_medium = classes.medium.model_minimal_medium(model, 'exchanges', growth_rate)
    else:
        minimal_medium = None

    # test for amino acid auxotrophies
    if test_aa_auxotrophies:
        auxo_results = auxotrophies_simulation(model, media_list, minimal_medium)

    # make and return report
    return classes.reports.GrowthAnalysisReport(growth_sim_results, minimal_medium, auxo_results)


# pathway analysis with KEGG
# --------------------------
def kegg_pathway_analysis(model):
    """Analyse the pathways that are covered by the model.

    The analysis is based on the KEGG pathway classification and the available
    KEGG pathway identifiers present in the model.
    Note: one reaction can have multiple pathway identifiers associated with it.

    :param model: The GEM to be analysed
    :type model: cobra.Model

    :returns: The results of the analysis.
    :rtype: PathwayAnalysisReport
    """

    report = classes.reports.PathwayAnalysisReport(total_reac=len(model.reactions))

    pathways = dict()
    counter = 0
    for r in model.reactions:
        if 'kegg.pathway' in r.annotation.keys():
            counter += 1
            anno = r.annotation['kegg.pathway']
            if isinstance(anno,str):
                anno = re.sub(r'^[a-z]*','',anno)
                if anno in pathways:
                    pathways[anno] += 1
                else:
                    pathways[anno] = 1
            else:
                for x in anno:
                    x = re.sub(r'^[a-z]*','',x)
                    if x in pathways:
                        pathways[x] += 1
                    else:
                        pathways[x] = 1

    report.kegg_count = counter

    global_map = {}
    over_map = {}
    rest = {}
    for k,v in pathways.items():
        if k.startswith('011'):
            global_map[k] = v
        elif k.startswith('012'):
            over_map[k] = v
        else:
            rest[k] = v

    report.kegg_global = global_map
    report.kegg_over = over_map
    report.kegg_paths = rest

    return report


# run this part
# -------------

def run(model_path, dir, pc_model_path=None, pc_based_on='id', db_path=None, load_media=['LB','M9'], change_to_aerobic=None,
        change_to_anaerobic=None, add_casamino=None, external_media=None,
        test_min_medium=False, growth_rate=0.1,
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
    model = util.io.read_model_cobra(model_path)

    # try to set objective to growth
    model.objective = util.cobra_models.find_growth_obj_func(model)

    # ------------------
    # general statistics
    # ------------------

    print('\n# ------------------\n# general statistics\n# ------------------')

    statistics_report = util.cobra_models.generate_statistics(model)
    statistics_report.save(F'{dir}05_analysis/')

    # -----------------
    # pan-core analysis
    # -----------------

    print('\n# ------------------\n# pan-core analysis\n# ------------------')

    if pc_model_path:
        pc_model = util.io.read_model_cobra(pc_model_path)
        pan_core_report = util.cobra_models.compare_to_core_pan(model, pc_model, pc_based_on)
        pan_core_report.save(F'{dir}05_analysis/')

    # ----------------
    # pathway analysis
    # ----------------

    print('\n# -----------------\n# pathway analysis\n# -----------------')

    if pathway:
        pathway_report = kegg_pathway_analysis(model)
        pathway_report.save(F'{dir}05_analysis/')

    # ---------------
    # growth analysis
    # ---------------

    print('\n# ---------------\n# growth analysis\n# ---------------')

    if load_media or external_media:
        growth_report = analyse_growth(model, db_path, load_media, change_to_aerobic,
                           change_to_anaerobic, add_casamino, external_media,
                           test_min_medium, growth_rate,
                           test_aa_auxotrophies)
        growth_report.save(F'{dir}05_analysis/')

    total_time_e = time.time()
    print(F'total runtime: {total_time_e-total_time_s}')
