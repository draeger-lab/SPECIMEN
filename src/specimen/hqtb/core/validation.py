"""Validate a model (part 4 of the workflow).

Implemented tests in include:
- cobra/sbml check using cobrapy
"""

__author__ = 'Carolin Brune'

################################################################################
# requirements
################################################################################

import cobra
import pprint
import time

from pathlib import Path
from typing import Literal

################################################################################
# functions
################################################################################

# @TODO: add more validation options.´

def run(dir:str, model_path:str, tests:None|Literal['cobra']=None, run_all:bool=True):
    """SPECIMEN Step 4: Validate the model.

    Args:
        - dir (str): 
            Path to the output directory.
        - model_path (str): 
            Path to the model to be validated
        - tests (None | Literal['cobra'], optional): 
            Tests to perform. 
            Defaults to None.
        - run_all (bool, optional): 
            Run al available tests. If True, overwrites
            the previous parameter.
            Defaults to True.
    """

    total_time_s = time.time()

    print('\nvalidation\n################################################################################\n')

    # -----------------------
    # create output directory
    # -----------------------

    try:
        Path(dir,"04_validation").mkdir(parents=True, exist_ok=False)
        print(F'Creating new directory {str(Path(dir,"04_validation"))}')
    except FileExistsError:
        print('Given directory already has required structure.')

    # --------------
    # validate model
    # --------------

    if run_all or (tests and 'cobra' in tests):
        print('\n# ---------------------------------\n# validate model - cobra validation\n# ---------------------------------')
        start = time.time()

        # validate using cobra
        cobra_report = cobra.io.validate_sbml_model(model_path)
        pprint.pprint(cobra_report)

        end = time.time()
        print(F'\ttime: {end - start}s')

    total_time_e = time.time()
    print(F'total runtime: {total_time_e-total_time_s}')
