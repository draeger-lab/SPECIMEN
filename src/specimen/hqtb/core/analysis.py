"""Analyse a model (step 5 of the workflow)."""

__author__ = "Carolin Brune"

################################################################################
# requirements
################################################################################

import logging
import time

from pathlib import Path
from typing import Literal

from refinegems.classes.medium import load_media
from refinegems.analysis import growth
from refinegems.utility.io import load_model
from refinegems.analysis.core_pan import compare_to_core_pan
from refinegems.curation.pathways import kegg_pathway_analysis
from refinegems.utility.util import test_biomass_presence
from refinegems.developement.decorators import suppress_warning

from ...classes.reports import SpecimenModelInfoReport

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

# run this part
# -------------

@suppress_warning("invalid character '*' found in formula")
def run(
    model_path: str,
    dir: str,
    media_path: str = None,
    namespace: Literal["BiGG"] = "BiGG",
    pc_model_path: str = None,
    pc_based_on: Literal["id"] = "id",
    test_aa_auxotrophies: bool = True,
    pathway: bool = True,
):
    """SPECIMEN Step 5: Analyse the generated model.

    Args:
        - model_path (str):
            Path to the model.
        - dir (str):
            Path to the output directory.
        - media_path (str, optional):
            Path to a media config file.
            Using this enables growth simulation.
            Defaults to None.
        - namespace (Literal['BiGG'], optional):
            Namespace to work on.
            Defaults to 'BiGG'.
        - pc_model_path (str, optional):
            Path to a core-pan model. Defaults to None.
        - pc_based_on (Literal['id'], optional):
            How to compare the model to the core-pan model.
            Defaults to 'id'.
        - test_aa_auxotrophies (bool, optional):
            Option to enable the amino acid
            auxotrophy simulation. Defaults to True.
        - pathway (bool, optional):
            Optional to enable KEGG pathway analysis.
            Defaults to True.
    """

    total_time_s = time.time()

    genlogger.info(
        "\nanalysis\n################################################################################\n"
    )

    # -----------------------
    # create output directory
    # -----------------------

    try:
        Path(dir, "05_analysis").mkdir(parents=True, exist_ok=False)
        genlogger.info(f'Creating new directory {Path(dir,"05_analysis")}')
    except FileExistsError:
        genlogger.info("Given directory already has required structure.")
        
    # set path for logging file
    Path(dir, "05_analysis", "analysis.log").unlink(missing_ok=True)
    handler = logging.handlers.RotatingFileHandler(
        str(Path(dir, "05_analysis", "analysis.log")),
        mode="w",
        # maxBytes=1000,
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
    # redirect matplotlib logging
    mpllogger = logging.getLogger("matplotlib")
    mpllogger.addHandler(handler)
    mpllogger.propagate = False
    # redirect refinegems logging
    rglogger = logging.getLogger("refinegems")
    rglogger.addHandler(handler)
    rglogger.propagate = False

    # load model
    model = load_model(str(model_path), "cobra")

    # ------------------
    # general statistics
    # ------------------

    logger.info("\n# ------------------\n# general statistics\n# ------------------")

    statistics_report = SpecimenModelInfoReport(model)
    statistics_report.save(Path(dir, "05_analysis"))

    # -----------------
    # pan-core analysis
    # -----------------

    if pc_model_path:
        logger.info("\n# ------------------\n# pan-core analysis\n# ------------------")
        pc_model = load_model(pc_model_path, "cobra")
        pan_core_report = compare_to_core_pan(model, pc_model, pc_based_on)
        pan_core_report.save(Path(dir, "05_analysis"))

    # ----------------
    # pathway analysis
    # ----------------

    if pathway:
        logger.info("\n# -----------------\n# pathway analysis\n# -----------------")
        pathway_report = kegg_pathway_analysis(model)
        pathway_report.save(Path(dir, "05_analysis"))

    # ---------------
    # growth analysis
    # ---------------

    if media_path:
        logger.info("\n# ---------------\n# growth analysis\n# ---------------")

        # try to set objective to growth
        growth_func_list = test_biomass_presence(model)
        if growth_func_list:
            # independently of how many growth functions are found, the first one will be used
            model.objective = growth_func_list[0]
            # simulate growth on different media
            growth_report = growth.growth_analysis(
                model, media_path, namespace=namespace, retrieve="report"
            )
            growth_report.save(Path(dir, "05_analysis"))

        else:
            logger.warning(
                "No growth/biomass function detected, growth simulation will be skipped."
            )

        # test auxotrophies
        if test_aa_auxotrophies:
            media_list = load_media(media_path)
            auxo_report = growth.test_auxotrophies(
                model, media_list[0], media_list[1], namespace
            )
            auxo_report.save(Path(dir, "05_analysis"))

    total_time_e = time.time()
    logger.info(f"total runtime: {total_time_e-total_time_s}")
