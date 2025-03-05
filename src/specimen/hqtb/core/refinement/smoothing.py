"""Perform part 4 of the refinement: smoothing.
"""

__author__ = 'Carolin Brune'

# @TODO logging 

################################################################################
# requirements
################################################################################

import time
import tempfile
import os
from pathlib import Path
from typing import Literal
import warnings

import cobra
from cobra import Reaction

from refinegems.curation.biomass import check_normalise_biomass
from refinegems.classes import egcs
from refinegems.utility.connections import adjust_BOF, perform_mcc, run_memote
from refinegems.utility.util import test_biomass_presence

################################################################################
# variables
################################################################################

################################################################################
# functions
################################################################################

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
    print('\n# ---------------------------------------------\n# # check for energy generating cycles\n# ---------------------------------------------')
    start = time.time()

    match egc_solver:
        # greedy solver
        case 'greedy':
            print('Using GreedyEGCSolver...')
            solver = egcs.GreedyEGCSolver()
            results = solver.solve_egcs(model,namespace=namespace) # @NOTE automatically uses c,p as compartments 
            if results:
                for k,v in results.items():
                    print(f'\t{k}: {v}')
        
        # no solver = EGCs will only be reported
        case _:
            solver = egcs.EGCSolver()
            print(f'\tFound EGCs:\n')
            print(f'\t{solver.find_egcs(model,with_reacs=True,namespace=namespace)}') # @NOTE automatically uses c,p as compartments 

    end = time.time()
    print(F'\ttime: {end - start}s')

    # --------------------------
    # biomass objective function
    # --------------------------
    # adjust the BOF to the current genome

    print('\n# ----------\n# adjust BOF\n# ----------')
    start = time.time()

    with tempfile.NamedTemporaryFile(suffix='.xml', delete=False) as temp_model:
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
            model.reactions.get_by_id(bof_reac).reaction = adjust_BOF(genome, temp_model.name, model, dna_weight_frac, ion_weight_frac)
       
            # optimise BOF(s)
            model = check_normalise_biomass(model)
    os.remove(temp_model.name)

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
        memote_path = str(Path(dir,'step4-smoothing',model_name+'.html'))
        run_memote(model, 'html', return_res=False, save_res=memote_path, verbose=True)

