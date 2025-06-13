"""Validate a model (step 4 of the workflow).

Implemented tests in include:
- cobra/sbml check using cobrapy
"""

__author__ = "Carolin Brune"

################################################################################
# requirements
################################################################################

import cobra
import logging
import pprint
import time
import warnings

import pandas as pd

from pathlib import Path
from typing import Union

from ...util.util import run_ModelPolisher

################################################################################
# setup logging
################################################################################
# general logging
genlogger = logging.getLogger(__name__)
# internal logger with logging file

logger = logging.getLogger(__name__ + "-intern")
logger.setLevel(logging.DEBUG)
logger.propagate = False

################################################################################
# functions
################################################################################


def run(
    dir: str,
    model_path: str,
    tests: Union[None, str, list] = None,
    run_all: bool = True,
):
    """SPECIMEN Step 4: Validate the model.
    
    Included tests (name : description):
    - modelpolisher: Semantic control and BiGG annotation fixing with ModelPolisher
    - cobra: SBML validation using COBRApy

    Args:
        - dir (str):
            Path to the output directory.
        - model_path (str):
            Path to the model to be validated
        - tests (Union[None, str, list], optional):
            Tests to perform.
            If the test name is either in a string or an element in a list, 
            the corresponding test will be run.
            Defaults to None.
        - run_all (bool, optional):
            Run al available tests. If True, overwrites
            the previous parameter.
            Defaults to True.
    """

    total_time_s = time.time()

    # -----------------------^
    # create output directory
    # -----------------------

    try:
        Path(dir, "04_validation").mkdir(parents=True, exist_ok=False)
        genlogger.info(f'Creating new directory {str(Path(dir,"04_validation"))}')
    except FileExistsError:
        genlogger.info("Given directory already has required structure.")

    # -----------------
    # fine tune logging
    # -----------------
    # interal logging
    Path(dir, "04_validation", "validation.log").unlink(missing_ok=True)
    handler = logging.handlers.RotatingFileHandler(
        str(Path(dir, "04_validation", "validation.log")),
        mode="w",
        backupCount=10,
        encoding="utf-8",
        delay=0,
    )
    handler.setFormatter(
        logging.Formatter(
            "{levelname} \t {name} \t {message}",
            style="{",
        )
    )
    logger.addHandler(handler)
    # redirect cobrapy logging
    cobralogger = logging.getLogger("cobra")
    cobralogger.addHandler(handler)
    cobralogger.propagate = False

    # --------------
    # validate model
    # --------------
    
    logger.info(
        "\nvalidation\n################################################################################\n"
    )
    
    # generalise input 
    match tests:
        case None:
            pass
        case str():
            tests = tests.lower()
        case list():
            tests = [t.lower() for t in tests]
        case _:
            warnings.warn(f"Tests parameter must be of type str or list, got {type(tests)}. Setting to None.")
            tests = None
            
    # ModelPolisher 
    # -------------
    
    if run_all or (tests and "modelpolisher" in tests):
        logger.info(
            "\n"
            "# -------------\n"
            "# ModelPolisher\n"
            "# -------------"
        )
        start = time.time()
        
        # generate specific directory for ModelPolisher output
        try:
            Path(dir, "04_validation", "modelpolisher").mkdir(parents=True, exist_ok=False)
            logger.info(f'Creating new directory {str(Path(dir,"04_validation", "modelpolisher"))}')
        except FileExistsError:
            logger.info("Given directory already has required structure.")
        
        # setting ModelPolisher params
        config_mp = {
            "allow-model-to-be-saved-on-server": False,
            "fixing": {"dont-fix": False},
            "annotation": {
                "bigg": {
                    "annotate-with-bigg": True,
                    "include-any-uri": False,
                }
            },
        }

        # running ModelPolisher
        result = run_ModelPolisher(model_path, config_mp)
        
        if result:
            # saving results files
            pd.DataFrame(result["diff"]).to_csv(
                Path(dir, "04_validation", "modelpolisher", "diff_mp.csv"),
                sep=";",
                header=False,
            )
            pd.DataFrame(result["pre_validation"]).to_csv(
                Path(dir, "04_validation", "modelpolisher", "pre_validation.csv"),
                sep=";",
                header=True,
            )
            pd.DataFrame(result["post_validation"]).to_csv(
                Path(dir, "04_validation", "modelpolisher", "post_validation.csv"),
                sep=";",
                header=True,
            )

            # save model
            # @ASK does this truly work? Test as soon as MP is up and running
            # From result the model object can be extracted like, I think: result["polished_document"].getModel()
            model_polisher_model_path = Path(dir, "04_validation", f"{Path(model_path).stem}_after_mp.xml")
            
            model_path = model_polisher_model_path
        
        end = time.time()
        logger.info(f"\ttime: {end - start}s")
    
    # COBRApy 
    # -------
    if run_all or (tests and "cobra" in tests):
        logger.info(
            "\n"
            "# ------------------\n"
            "# COBRApy validation\n"
            "# ------------------"
        )
        start = time.time()

        # validate using cobra
        cobra_report = cobra.io.validate_sbml_model(model_path)
        with open(
            Path(dir, "04_validation", "cobrapy-validation.txt"), "w"
        ) as cpyval_file:
            pprint.pprint(cobra_report, stream=cpyval_file)

        end = time.time()
        logger.info(f"\ttime: {end - start}s")

    total_time_e = time.time()
    logger.info(f"total runtime: {total_time_e-total_time_s}")

    # restore cobrapy logging behaviour
    cobralogger.handlers.clear()
    cobralogger.propagate = False
