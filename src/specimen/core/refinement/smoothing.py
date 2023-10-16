"""Perform part 4 of the refinement: smoothing.
"""

__author__ = 'Carolin Brune'

################################################################################
# requirements
################################################################################

import subprocess
import time
import tempfile
from pathlib import Path
import warnings

from BOFdat import step1
from BOFdat import step2
from BOFdat.util import update
from BOFdat.util.update import determine_coefficients
from MCC import MassChargeCuration

import pandas as pd
import cobra

from ... import util

# further required programs:
#        - MEMOTE,  tested with version 0.13.0
#        - BOFdat
#        - MassChargeCuration

# note:
#    for BOFdat to run correctly, you need to change 'solution.f' to 'solution.objective_value'
#    in the coenzymes_and_ions.py file of BOFdat

################################################################################
# variables
################################################################################

dissipation_reactions = {'DR_ATP':{'atp_c':-1.0,'h2o_c':-1.0,'h_c':1.0,'adp_c':1.0,'pi_c':1.0},
                        'DR_CTP':{'ctp_c':-1.0,'h2o_c':-1.0,'h_c':1.0,'cdp_c':1.0,'pi_c':1.0},
                        'DR_GTP':{'gtp_c':-1.0,'h2o_c':-1.0,'h_c':1.0,'gdp_c':1.0,'pi_c':1.0},
                        'DR_UTP':{'utp_c':-1.0,'h2o_c':-1.0,'h_c':1.0,'udp_c':1.0,'pi_c':1.0},
                        'DR_ITP':{'itp_c':-1.0,'h2o_c':-1.0,'h_c':1.0,'idp_c':1.0,'pi_c':1.0},
                        'DR_NADH':{'nadh_c':-1.0,'h_c':1.0,'nad_c':1.0},
                        'DR_NADPH':{'nadph_c':-1.0,'h_c':1.0,'nadp_c':1.0},
                        'DR_FADH2':{'fadh2_c':-1.0,'h_c':2.0,'fad_c':1.0},
                        'DR_FMNH2':{'fmnh2_c':-1.0,'h_c':2.0,'fmn_c':1.0},
                        'DR_Q8H2':{'q8h2_c':-1.0,'h_c':2.0,'q8_c':1.0},
                        'DR_MQL8':{'mql8_c':-1.0,'h_c':2.0,'mqn8_c':1.0},
                        'DR_DMMQL8':{'2dmmql8_c':-1.0,'h_c':2.0,'2dmmq8_c':1.0},
                        'DR_ACCOA':{'h2o_c':-1.0,'accoa_c':-1.0,'h_c':1.0,'ac_c':1.0,'coa_c':1.0},
                        'DR_GLU':{'h2o_c':-1.0,'glu__L_c':-1.0,'h_c':2.0,'akg_c':1.0,'nh4_c':1.0},
                        'DR_PROTON':{'h_p':-1.0,'h_c':1.0}}

################################################################################
# functions
################################################################################

def perform_mcc(model, dir, apply = True):
    """Use the MassChargeCuration tool on the model. Either apply it directly
    or on a copy.

    :param model: The model to use the tool on.
    :type  model: cobra.Model
    :param dir:   Path of the directory to save MCC output in.
    :type  dir:   string
    :param apply: Option, if MCC should be used directly on the model or not.
    :type  apply: bool, true if MCC should be used directly on the model.
    :returns:     The model (updated or not).
    :rtype:       cobra.Model
    """

    # make temporary directory to save files for MCC in
    with tempfile.TemporaryDirectory() as temp:

        # use MCC
        if apply:
            # update model
            balancer = MassChargeCuration(model, update_ids = False, data_path=temp)
        else:
            # do not change original model
            model_copy = model.copy()
            balancer = MassChargeCuration(model_copy, update_ids = False, data_path=temp)

    # save reports
    balancer.generate_reaction_report(F'{dir}{model.id}_mcc_reactions')
    balancer.generate_metabolite_report(F'{dir}{model.id}_mcc_metabolites')
    balancer.generate_visual_report(F'{dir}{model.id}_mcc_visual')

    return model


def check_for_cycles(model):
    """Check a given metabolic model for futile cycles and energy generating cycles.
    As the correction of those is highly dependant on biological evidence, the correction
    has to be performed manually after using this function to check the existence and possible
    reasons for cycles.

    :param model: The genome–scale metabolic model.
    :type  model: cobra.Model
    """

    external_compartment = cobra.medium.find_external_compartment(model)
    print(F'\tcompartment found as external compartment: {external_compartment}')
    for dis_reac_id, dis_reac_dict in dissipation_reactions.items():

        # make sure not to wirk on the model directly while checking for cycles
        with model as test_model:

            # remove medium
            test_model.medium = {}

            # create dissipation reaction
            dis_reac = cobra.Reaction(dis_reac_id, name = F'Dissipation reaction {dis_reac_id}')

            # find metabolites in model
            metabolites = {}
            missing_metabolite = False
            for m,f in dis_reac_dict.items():
                try:
                    metabolites[test_model.metabolites.get_by_id(m)] = f
                except:
                    # if metabolites not found ...
                    #@ TODO
                    missing_metabolite = True
            # check if all metabolites have been found in the model
            if missing_metabolite:
                print(F'\tMissing metabolite in {dis_reac_id}, reaction skipped.')
                continue
            # add metabolites to the reaction
            dis_reac.add_metabolites(metabolites)
            # add reaction to model
            test_model.add_reactions([dis_reac])

            # set reaction as model objective function
            test_model.objective = dis_reac.id

            # check for futile cycles
            solution_futile = test_model.optimize().objective_value
            if solution_futile != 0.0:
                print(F'\tFound futile cycle for {dis_reac.id} with {solution_futile}')
                for reac in test_model.reactions:
                    pass
                # --------
                # problem
                #@TODO
                # print suggestion on what to do manually
                # or any good idea to somehow influence it
                # --------

            # reassign boundaries depending on reaction type
            for reac in test_model.reactions:
                # exchange reactions
                if cobra.medium.is_boundary_type(reac,'exchange',external_compartment):
                    reac.upper_bound = 0.0
                    reac.lower_bound = 0.0
                # reversible reactions
                elif reac.reversibility:
                    reac.upper_bound = 1.0
                    reac.lower_bound = -1.0
                # irreversible reactions (remaining)
                elif not reac.reversibility:
                    reac.upper_bound = 1.0
                    reac.lower_bound = 0.0
                else:
                    pass

            # reset boundaries for dissipation reaction
            test_model.reactions.get_by_id(dis_reac.id).upper_bound = 1000

            # check for EGC (energy generating cycles)
            test_model.optimize()
            for reac in test_model.reactions:
                if reac.flux != 0.0:
                    print(F'\tFound EGC: flux = {reac.flux} for {dis_reac.id} at {reac.id}')
                    # --------
                    # problem
                    #@TODO
                    # print suggestion on what to do manually
                    # or any good idea to somehow influence it
                    # --------


def adjust_BOF(genome, model_file, model, dna_weight_fraction, weight_frac):
    """Adjust the model's BOF using BOFdat. Currently implemented are step 1
    dna coefficients and step 2.

    :param genome:              Path to the genome (e.g. .fna) FASTA file.
    :type  genome:              string
    :param model_file:          Path to the sbml (.xml) file of the model.
    :param model_file:          string
    :param model:               The loaded genome-scale metabolic model (from the string above).
    :type  model:               cobra.Model
    :param dna_weight_fraction: DNA weight fraction for BOF step 1.
    :type  dna_weight_fraction: float
    :param weight_frac:         Weight fraction for the second step of BOFdat (enzymes and ions)
    :type  weight_frac:         float
    :returns:                   The updated BOF reaction.
    :rtype:                     string
    """

    # BOFdat step 1:
    # --------------
    # dna coefficients
    dna_coefficients = step1.generate_dna_coefficients(genome,model_file,DNA_WEIGHT_FRACTION=dna_weight_fraction)
    bd_step1 = {}
    for m in dna_coefficients:
        bd_step1[m.id] = dna_coefficients[m]

    # ...........................
    # @TODO
    #    if time permits or needed, options for more coefficients can be added
    # ...........................

    # BOFdat step 2:
    # --------------
    # find inorganic ions
    selected_metabolites = step2.find_coenzymes_and_ions(model_file)
    # determine coefficients
    bd_step2 = determine_coefficients(selected_metabolites,model,weight_frac)
    bd_step2.update(bd_step1)

    # update BOF
    # ----------
    #  retrieve previously used BOF
    growth_objname = util.cobra_models.find_growth_obj_func(model)
    objective_list = model.reactions.get_by_id(growth_objname).reaction.split(' ')
    # ...............................................................
    objective_reactant = {}
    objective_product = {}
    product = False
    # get reactants, product and factors from equation
    for s in objective_list:
        if s == '+':
            continue
        elif '>' in s:
            product = True
        elif '.' in s:
            factor = s
        else:
            metabolite = s
            if product:
                objective_product[s] = factor
            else:
                objective_reactant[s] = factor

    # update BOF information with data from BOFdat
    for m in bd_step2:
        if bd_step2[m] < 0:
            objective_reactant[m] = bd_step2[m]*(-1)
        else:
            objective_product[m] = bd_step2[m]

    # check if values are within range
    #     (sum(reactant * reactant_mol_weight) - sum(product * mol_weight)) / 1000 = x
    #     acceptance if solution  1-x < -1e-06 or 1-x > 1e-03
    t = 0
    for m,x in objective_reactant.items():
        t += float(x) * float(model.metabolites.get_by_id(m).formula_weight)/1000
    for m,x in objective_product.items():
        t -= float(x) * float(model.metabolites.get_by_id(m).formula_weight)/1000

    i = 0
    # if not, normalise them until they are
    # or 10 rounds of correction have been run
    while 1-t < -1e-06 or 1-t > 1e-03 or i >= 10:
        y = t
        t = 0
        for m,x in objective_reactant.items():
            objective_reactant[m] = float(x)/y
            t += (float(x)/y) * float(model.metabolites.get_by_id(m).formula_weight)/1000
        for m,x in objective_product.items():
            objective_product[m] = float(x)/y
            t -= (float(x)/y) * float(model.metabolites.get_by_id(m).formula_weight)/1000
        i += 1

    # create new objective function
    new_objective = ' + '.join('{} {}'.format(value, key) for key, value in objective_reactant.items())
    new_objective += ' --> '
    new_objective += ' + '.join('{} {}'.format(value, key) for key, value in objective_product.items())

    return new_objective


def run(genome,model,dir,mcc='skip',dna_weight_frac=0.023,ion_weight_frac=0.05, memote=False):
    """Perform the fourth step of the refinement, smoothing, on a model.

    The fourth step of the refinment, smoothing, includes:
    - Mass-Charge curation
    - checking for futile and energy generating cycles
    - adjusting Biomass objective function parameters based on genome

    :param genome: Path to the genome FASTA (e.g. .fna) file of your genome.
    :type genome: string
    :param model: Path to the model. Should be sbml.
    :type model: string
    :param dir: Path of the output directory.
    :type dir: string

    :param mcc: Option for the Mass-Charge curation.
        Can be 'skip', 'apply' (directly use MCC on model) or 'extra' (only generate information).
        Default is 'skip'.
    :type mcc: string
    :param dna_weight_frac: DNA weight fraction to use for BOFdat.
        Default is 0.023 for Klebsiella pneumoniae based on Liao et al.
    :type dna_weight_frac: float
    :param ion_weight_frac: Ion weight fraction to use for BOFdat.
        Default is the BOFdat default.
    :type ion_weight_frac: float
    :param memote: Option to run memote after the refinement.
        Default is False.
    :type memote: bool
    """


    print('\nrefinement step 4: smoothing\n################################################################################\n')

    # -----------------------
    # create output directory
    # -----------------------

    # make sure given directory path ends with '/'
    if not dir.endswith('/'):
        dir = dir + '/'

    try:
        Path(F"{dir}step4-smoothing/").mkdir(parents=True, exist_ok=False)
        print(F'Creating new directory {F"{dir}step4-smoothing/"}')
    except FileExistsError:
        print('Given directory already has required structure.')

    try:
        Path(F"{dir}manual_curation/").mkdir(parents=True, exist_ok=False)
        print(F'Creating new directory {F"{dir}manual_curation/"}')
    except FileExistsError:
        print('Given directory already has required structure.')

    # ---------
    # load data
    # ---------
    model = cobra.io.read_sbml_model(model)

    # ---------------
    # mass and charge
    # ---------------

    if mcc == 'apply':
        print('\n# ----------------------------------\n# mass and charge curation (applied)\n# ----------------------------------')
        start = time.time()
        model = perform_mcc(model,F"{dir}step4-smoothing/")
        end = time.time()
        print(F'\ttime: {end - start}s')
    elif mcc == 'extra':
        print('\n# --------------------------------\n# mass and charge curation (extra)\n# --------------------------------')
        start = time.time()
        model = perform_mcc(model,F"{dir}manual_curation/",False)
        end = time.time()
        print(F'\ttime: {end - start}s')
    elif mcc == 'skip':
        print('\n# ------------------------\n# mass and charge curation\n# ------------------------\n\tskipped')
    else:
        warnings.warn(F'Unknown option {mcc} for Mass and Charge Curation. Usage of MCC will be skipped.')

    # ---------------------------------------------
    # check for futile and energy generating cycles
    # ---------------------------------------------
    # check if some exist
    # not implemented: solve them, as it is heavily dependant on biology (hard to perform automated)
    print('\n# ---------------------------------------------\n# # check for futile and energy generating cycles\n# ---------------------------------------------')
    start = time.time()

    check_for_cycles(model)

    end = time.time()
    print(F'\ttime: {end - start}s')

    # --------------------------
    # biomass objective function
    # --------------------------
    # adjust the BOF to the current genome

    print('\n# ----------\n# adjust BOF\n# ----------')
    start = time.time()

    with tempfile.NamedTemporaryFile(suffix='.xml') as temp_model:
        # generate an up-to-date model xml-file
        cobra.io.write_sbml_model(model,temp_model.name)
        # update BOF
        model.reactions.get_by_id(util.cobra_models.find_growth_obj_func(model)).reaction = adjust_BOF(genome, temp_model.name, model, dna_weight_frac, ion_weight_frac)

    end = time.time()
    print(F'\ttime: {end - start}s')

    # ----------------
    # save final model
    # ----------------
    print('\n# ----------\n# save model\n# ----------')
    model_name = F'{model.id}_smooth'
    outname = F"{dir}step4-smoothing/{model_name}.xml"
    print(F'\tsaving to: {outname}')
    cobra.io.write_sbml_model(model,outname)

    # ---------------------------------
    # assess model quality using memote
    # ---------------------------------

    if memote:
        print('\n# -------------------\n# analyse with MEMOTE\n# -------------------')
        start = time.time()
        draft_path = outname.replace(" ", "\ ")
        memote_path = F'{dir}step4-smoothing/{model_name}.html'.replace(" ", "\ ")
        subprocess.run([F'memote report snapshot --filename {memote_path} {draft_path}'], shell=True)
        end = time.time()
        print(F'\ttotal time: {end - start}s')
