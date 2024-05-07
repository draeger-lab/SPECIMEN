"""Perform part 4 of the refinement: smoothing.
"""

__author__ = 'Carolin Brune'

################################################################################
# requirements
################################################################################

import time
import tempfile
from pathlib import Path
from typing import Literal
import warnings

from BOFdat import step1
from BOFdat import step2
from BOFdat.util import update
from BOFdat.util.update import determine_coefficients
from MCC import MassChargeCuration

import cobra
from cobra import Reaction

from refinegems.curation.biomass import test_biomass_presence, check_normalise_biomass
from refinegems.analysis.investigate import run_memote
from refinegems.classes import egcs

# further required programs:
#        - BOFdat
#        - MassChargeCuration

# note:
#    for BOFdat to run correctly, you need to change 'solution.f' to 'solution.objective_value'
#    in the coenzymes_and_ions.py file of BOFdat

################################################################################
# variables
################################################################################

################################################################################
# functions
################################################################################

def perform_mcc(model: cobra.Model, dir: str, apply:bool = True) -> cobra.Model:
    """Run the MassChargeCuration toll on the model and optionally directly apply 
    the solution.

    Args:
        - model (cobra.Model): 
            The model to use the tool on.
        - dir (str): 
            Path of the directory to save MCC output in.
        - apply (bool, optional): 
            If True, model is directly updated with the results. 
            Defaults to True.

    Returns:
        cobra.Model: 
            The model (updated or not)
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
    balancer.generate_reaction_report(Path(dir,model.id+'_mcc_reactions'))
    balancer.generate_metabolite_report(Path(dir,model.id+'_mcc_metabolites'))
    balancer.generate_visual_report(Path(dir,model.id+'_mcc_visual'))

    return model


def adjust_BOF(genome:str, model_file:str, model:cobra.Model, dna_weight_fraction:float, weight_frac:float) -> str:
    """Adjust the model's BOF using BOFdat. Currently implemented are step 1
    DNA coefficients and step 2.

    Args:
        - genome (str): 
            Path to the genome (e.g. .fna) FASTA file.
        - model_file (str): 
            Path to the sbml (.xml) file of the model.
        - model (cobra.Model): 
            The genome-scale metabolic model (from the string above), loaded with COBRApy.
        - dna_weight_fraction (float): 
            DNA weight fraction for BOF step 1.
        - weight_frac (float): 
            Weight fraction for the second step of BOFdat (enzymes and ions)

    Returns:
        str: 
            The updated BOF reaction as a reaction string.
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
    growth_func_list = test_biomass_presence(model)
    if len(growth_func_list) == 1: 
        objective_list = model.reactions.get_by_id(growth_func_list[0]).reaction.split(' ')
    elif len(growth_func_list) > 1:
        mes = f'Multiple BOFs found. Using {growth_func_list[0]} for BOF adjustment.'
        warnings.warn(mes,category=UserWarning)
        objective_list = model.reactions.get_by_id(growth_func_list[0]).reaction.split(' ')
    # else not needed, as if there is no BOF, a new one will be created

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

    # create new objective function
    new_objective = ' + '.join('{} {}'.format(value, key) for key, value in objective_reactant.items())
    new_objective += ' --> '
    new_objective += ' + '.join('{} {}'.format(value, key) for key, value in objective_product.items())

    return new_objective


def run(genome:str,model:str,dir:str,mcc='skip',
        egc_solver:None|Literal['greedy']=None,
        namespace:Literal['BiGG']='BiGG',
        dna_weight_frac=0.023,ion_weight_frac=0.05, 
        memote=False):
    """Perform the fourth step of the refinement, smoothing, on a model.

    The fourth step of the refinment, smoothing, includes:

    - Mass-Charge curation
    - checking energy generating cycles
    - adjusting Biomass objective function parameters based on genome

    Args:
        - genome (str): 
            Path to the genome FASTA (e.g. .fna) file of your genome.
        - model (str): 
            Path to the model. Should be sbml-format.
        - dir (str): 
            Path of the output directory.
        - mcc (str, optional): 
            Option for the Mass-Charge curation.
            Can be 'skip', 'apply' (directly use MCC on model) or 'extra' (only generate information). 
            Defaults to 'skip'.
        - egc_solver (None | Literal['greedy'], optional): 
            If given, uses the option to solve EGCs.
            Current options include greedy (Greedy Solver).
            Defaults to 'greedy'.
        - namespace (Literal['BiGG'], optional): 
            Namespace to use for the model.
            Defaults to 'BiGG'.
        - dna_weight_frac (float, optional): 
            DNA weight fraction to use for BOFdat.
            Default is 0.023 for Klebsiella pneumoniae based on Liao et al.
        - ion_weight_frac (float, optional): 
            Ion weight fraction to use for BOFdat.
            Defaults to 0.05.
        - memote (bool, optional): 
            Option to run memote after the refinement.
            Defaults to False.
    """
    
    print('\nrefinement step 4: smoothing\n################################################################################\n')

    # -----------------------
    # create output directory
    # -----------------------

    try:
        Path(dir,"step4-smoothing").mkdir(parents=True, exist_ok=False)
        print(F'Creating new directory {str(Path(dir,"step4-smoothing"))}')
    except FileExistsError:
        print('Given directory already has required structure.')

    try:
        Path(dir,"manual_curation").mkdir(parents=True, exist_ok=False)
        print(F'Creating new directory {str(Path(dir,"manual_curation"))}')
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
        model = perform_mcc(model,Path(dir,"step4-smoothing"))
        end = time.time()
        print(F'\ttime: {end - start}s')
    elif mcc == 'extra':
        print('\n# --------------------------------\n# mass and charge curation (extra)\n# --------------------------------')
        start = time.time()
        model = perform_mcc(model,Path(dir,"manual_curation"),False)
        end = time.time()
        print(F'\ttime: {end - start}s')
    elif mcc == 'skip':
        print('\n# ------------------------\n# mass and charge curation\n# ------------------------\n\tskipped')
    else:
        warnings.warn(F'Unknown option {mcc} for Mass and Charge Curation. Usage of MCC will be skipped.')

    # ----------------------------------
    # check for energy generating cycles
    # ----------------------------------
    # check if some exist
    # not implemented: solve them, as it is heavily dependant on biology (hard to perform automated)
    print('\n# ---------------------------------------------\n# # check for energy generating cycles\n# ---------------------------------------------')
    start = time.time()

    match egc_solver:
        # greedy solver
        case 'greedy':
            print('GreedyEGCSolver:')
            solver = egcs.GreedyEGCSolver()
            results = solver.solve_egcs(model,namespace=namespace) # @NOTE automatically uses c,p as compartments - maybe change later
            if results:
                for k,v in results.items():
                    print(f'\t{k}: {v}')
        
        # no solver = EGCs will only be reported
        case _:
            solver = egcs.EGCSolver()
            print(f'\tFound EGCs:\n')
            print(f'\t{solver.find_egcs()}')

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
        pos_bofs = test_biomass_presence(model)
        if pos_bofs:
            model.reactions.get_by_id(pos_bofs[0]).reaction = adjust_BOF(genome, temp_model.name, model, dna_weight_frac, ion_weight_frac)
            # optimise BOF(s)
            model = check_normalise_biomass(model)
        else:
            # create new BOF
            bof_reac = Reaction('Biomass_BOFdat')
            bof_reac.name = 'Biomass objective function created by BOFdat'
            model.add_reactions([bof_reac])
            model.reactions.get_by_id(pos_bofs[0]).reaction = adjust_BOF(genome, temp_model.name, model, dna_weight_frac, ion_weight_frac)
       
            # optimise BOF(s)
            model = check_normalise_biomass(model)


    end = time.time()
    print(F'\ttime: {end - start}s')

    # ----------------
    # save final model
    # ----------------
    print('\n# ----------\n# save model\n# ----------')
    model_name = F'{model.id}_smooth'
    outname = Path(dir,'step4-smoothing',model_name+".xml")
    print(F'\tsaving to: {outname}')
    cobra.io.write_sbml_model(model,outname)

    # ---------------------------------
    # assess model quality using memote
    # ---------------------------------

    if memote:
        memote_path = Path(dir,'step4-smoothing',model_name+'.html')
        run_memote(model, 'html', return_res=False, save_res=memote_path, verbose=True)

