"""Validate a model (part 4 of the workflow).

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

from pathlib import Path
from typing import Literal

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
    tests: None | Literal["cobra"] = None,
    run_all: bool = True,
):
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

    # -----------------------
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

    if run_all or (tests and "cobra" in tests):
        logger.info(
            "\n# ---------------------------------\n# validate model - cobra validation\n# ---------------------------------"
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
