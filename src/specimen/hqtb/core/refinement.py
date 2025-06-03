"""Function for the refinement step (step 3) of the HQTB workflow.

Refinement includes four main parts:
- Part 1: extension
- Part 2: cleanup
- Part 3: annotation
- Part 4: smoothing
"""

__author__ = "Carolin Brune"

################################################################################
# requirements
################################################################################

import cobra
import logging
import os
import time
import tempfile
import warnings

from cobra import Reaction
from libsbml import *
from pathlib import Path
from tempfile import NamedTemporaryFile
from tqdm import tqdm
from tqdm.auto import tqdm

tqdm.pandas()
from typing import Literal, Union

# refinegems
from refinegems.analysis.growth import load_media
from refinegems.classes import egcs
from refinegems.classes.gapfill import multiple_cobra_gapfill, GeneGapFiller
from refinegems.curation.curate import resolve_duplicates, check_direction
from refinegems.curation.biomass import check_normalise_biomass
from refinegems.curation.pathways import set_kegg_pathways
from refinegems.curation.miriam import polish_annotations
from refinegems.utility.connections import (
    adjust_BOF,
    perform_mcc,
    run_memote,
    run_SBOannotator
)
from refinegems.utility.io import load_model, write_model_to_file, parse_gff_for_cds, convert_cobra_to_libsbml
from refinegems.utility.util import test_biomass_presence

# further required programs:
#        - DIAMOND, tested with version 0.9.14 (works only for certain sensitivity mode)
#                   tested with version 2.1.8 (works for all sensitivity modes for that version)

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
# variables
################################################################################

GGF_REQS = {
    "prefix",
    "type_db",
    "fasta",
    "dmnd_db",
    "map_db",
    "mail",
    "check_NCBI",
    "threshold_add_reacs",
    "sens",
    "cov",
    "t",
    "pid",
    "formula_check",
    "exclude_dna",
    "exclude_rna",
    "gff",
} #: :meta:

################################################################################
# functions
################################################################################

# Part 1: extension
# -----------------

def extend(
    draft: str,
    gff: str,
    fasta: str,
    db: str,
    dir: str,
    # mapping to NCBI
    ncbi_mapping: Union[Path, str, None] = None,
    email: Union[None, str] = None,
    # params for DIAMOND
    sensitivity: Literal[
        "sensitive", "more-sensitive", "very-sensitive", "ultra-sensitive"
    ] = "more-sensitive",
    coverage: float = 95.0,
    pid: float = 90.0,
    threads: int = 2,
    # param for adding entities to model
    threshold_add_reacs: int = 5,
    prefix: str = "specimen",
    namespace: str = "BiGG",
    formula_check: Literal["none", "existence", "wildcard", "strict"] = "existence",
    exclude_dna: bool = True,
    exclude_rna: bool = True,
    # validation
    memote: bool = False,
):
    """Extend a draft model.

    Using a multitude of information from MetaNetX, KEGG, NCBI and more,
    extend a draft model by yet unadded genes ("missing genes") by mapping them against
    a database using DIAMOND.

    Args:
        - draft (str):
            Path to the draft model.
        - gff (str):
            GFF file of the organism of interest.
        - fasta (str):
            Path to the (protein) FASTA file containing the CDS sequences.
            
            .. warning:: 
                This FASTA needs to be in the extended GenBank format.  
                This can be downloaded from the NCBI FTB servers 
                (name usually contains `_tranlated_CDS`) or use 
                `mimic_genbank` from refineGEMs to create a similar format.
            
        - db (str):
            Path to the database used for running DIAMOND.
        - dir (str):
            Path to the directory for the output (directories).

        - ncbi_mapping (str):
            Path to the NCBI information mapping file. Optional, but recommended.
            Significantly decreases the runtime.
        - email (str):
            User's mail address for accessing information from NCBI.
            Will only be used, if ncbi_mapping is not set.
            Significantly increases the runtime.

        - sensitivity (Literal['sensitive', 'more-sensitive', 'very-sensitive','ultra-sensitive'], optional):
            Sensitivity mode for DIAMOND blastp run.
            Can be sensitive, more-sensitive, very-sensitive or ultra-sensitive.
            Defaults to 'more-sensitive'.
        - coverage (float, optional):
            Threshold value for the query coverage for DIAMOND.
            Defaults to 95.0.
        - pid (float, optional):
            PID (percentage identity value) to filter the blast hist by.
            Only hits equal or above the given value are kept.
            Defaults to 90.0.
        - threads (int, optional):
            Number of threads to be used for DIAMOND.
            Defaults to 2.

        - threshold_add_reacs (int, optional):
            Maximal number of reactions IDs allowed for one EC number mapping.
            If more mappings are found, the corresponding reeactions will not be
            added to the model due to low specifity.
            Defaults to 5.
        - prefix (str, optional):
            If the found attributes for a model entity cannot be used for
            building an ID, this prefix will be added to it.
            Defaults to 'specimen'
        - namespace (str, optional):
            Namespace to use for the added entities.
            Defaults to 'BiGG'.

            .. note::
                Based on the current developement status of refinegems,
                it is advised to

        - formula_check (Literal['none','existence','wildcard','strict'], optional)
            Level of chemical formula to be accespted for adding reactions.
            For more information, refer to the refineGems documentation.
            Defaults to 'existence'.

        - exclude_dna (bool, optional):
            Exclude reactions with DNA in their name when added.
            Defaults to True.
        - exclude_rna (bool, optional):
            Exclude reactions with RNA in their name when added.
            Defaults to True.

        - memote (bool, optional):
            Use memote on the extended model.
            Defaults to False.

    Raises:
        - ValueError: Unknown sensitivity mode.
    """

    if not sensitivity in [
        "sensitive",
        "more-sensitive",
        "very-sensitive",
        "ultra-sensitive",
    ]:
        raise ValueError(
            f"Unknown sensitivity mode {sensitivity}. Choose from: sensitive, more-sensitive, very-sensitive, ultra-sensitive"
        )

    logger.info(
        "\nrefinement step 1: extension\n################################################################################\n"
    )

    # -----------------------
    # create output directory
    # -----------------------

    try:
        Path(dir, "step1-extension").mkdir(parents=True, exist_ok=False)
        logger.info(f'Creating new directory {str(Path(dir,"step1-extension"))}')
    except FileExistsError:
        logger.info("Given directory already has required structure.")
        
    # set path for logging file
    Path(dir, "step1-extension", "extension.log").unlink(missing_ok=True)
    handler = logging.handlers.RotatingFileHandler(
        str(Path(dir, "step1-extension", "extension.log")),
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
    
    # -------------
    # add locus tag
    # -------------
    
    # get ncbiprotein + locus tag mapping
    gfftable = parse_gff_for_cds(
        gff, keep_attributes = {"locus_tag": "locus_tag", "protein_id": "ncbiprotein"}
    )
    
    # load model 
    draft_cobra = load_model(draft, "cobra")
    
    # identify, where locus tags are stored / can be stored
    add_label_via = None
    if len(set(g.id for g in draft_cobra.genes).union(set(gfftable['locus_tag']))) == 0:
        # case 1: ncbiprotein as ID
        # add locus tags as notes
        if 'ncbiprotein' in gfftable.columns:
            gfftable = gfftable.explode('ncbiprotein')
            add_label_via = 'notes'
            for g in draft_cobra.genes:
                hit = list(gfftable[gfftable['ncbiprotein']==g.id]['locus_tag'])
                if hit:
                    g.notes['locus_tag'] =  hit[0]
        else:
            raise KeyError('No NCBIprotein column detected. Possibly an error GFF file or draft model IDs. Please re-check your input.')

    else:
        # case 2: locus tag as ID
        add_label_via = 'id'
    
    # also generate libsbml instance
    # with locus tag as labels
    draft_libsbml = convert_cobra_to_libsbml(draft_cobra, add_label_via)
    
    # ----------------------
    # identify missing genes
    # ----------------------

    name = f"{draft_cobra.id}_extended"

    # set up GapFiller
    gp = GeneGapFiller()

    # identify missing genes
    gp.find_missing_genes(gff, draft_libsbml)

    if hasattr(gp, "missing_genes") and len(gp.missing_genes) >= 1:

        # --------------------------
        # identify missing reactions
        # --------------------------

        kwargs = {
            "sens": sensitivity,
            "cov": coverage,
            "t": threads,
            "pid": pid,
            "outdir": dir,
        }
        gp.find_missing_reactions(
            draft_cobra,
            prefix,
            "user",
            fasta,
            db,
            ncbi_mapping,
            email,
            False,
            threshold_add_reacs,
            **kwargs,
        )

        # ------------
        # extend model
        # ------------

        kwargs = {
            "namespace": namespace,
            "idprefix": prefix,
            "formula_check": formula_check,
            "exclude_dna": exclude_dna,
            "exclude_rna": exclude_rna,
        }
        extended_model = gp.fill_model(draft_libsbml, **kwargs)

        # save GapFiller report
        gp.report(Path(dir, "step1-extension"))

        # save model
        write_model_to_file(extended_model, Path(dir, "step1-extension", name + ".xml"))

        # ---------------------------------
        # assess model quality using memote
        # ---------------------------------

        if memote:
            logger.info("\nRunning memote ...\n------------------\n")
            memote_path = str(Path(dir, "step1-extension", name + ".html"))
            run_memote(
                draft, "html", return_res=False, save_res=memote_path, verbose=True
            )
    else:
        logger.warning(
            "\tNo missing genes found. The first refinement step, extension, will be skipped. The model will still be re-saved under the new name."
        )
        # save model
        write_model_to_file(draft_libsbml, Path(dir, "step1-extension", name + ".xml"))


# Part 2: cleanup
# ---------------

def cleanup(
    model: str,
    dir: str,
    biocyc_db: str = None,
    run_gene_gapfiller: Union[None, dict] = None,
    check_dupl_reac: bool = False,
    check_dupl_meta: bool = "default",
    remove_unused_meta: bool = False,
    remove_dupl_reac: bool = False,
    remove_dupl_meta: bool = False,
    universal: str = None,
    media_path: str = None,
    namespace: Literal["BiGG"] = "BiGG",
    growth_threshold: float = 0.05,
    iterations: int = 3,
    chunk_size: int = 10000,
    memote: bool = False,
):
    """Perform the second refinement step, cleanup, on a model.

    The second refinement step resolves the following issues:

    (1) (optional) checking direction of reactions with BioCyc
    (2) (optional) gapfilling
        (a) using refineGEMs GeneGapFiller
        (b) using cobra + universal model with reactions + a set of media
    (3) find and/or resolve duplicates (reactions and metabolites)

    Args:
        - model (str):
            The Path to an sbml model.
        - dir (str):
            Path to the directory of the output.
        - biocyc_db (str, optional):
            Path to the BioCyc/MetaCyc reaction database file.
            Defaults to None, which leads to skipping the direction check.
        - run_gene_gapfiller (Union[None,dict], optional):
            If a dictionary is given, tries to run the GeneGapFiller.
            If set to None, the gap-filling will be skipped
            The dictionary needs to contain the keys saved in :py:data:`~specimen.hqtb.core.refinement.GGF_REQS`.
            Defaults to None.
        - check_dupl_reac (bool, optional):
            Option to check for duplicate reactions.
            Defaults to False.
        - check_dupl_meta (bool, optional):
            Option to check for duplicate metabolites.
            Defaults to 'default', which checks based on MetaNetX first.
            Further options include 'skip' and 'exhaustive' (check for all possibilities).
        - remove_unused_meta (bool, optional):
            Option to remove unused metabolites from the model.
            Defaults to False.
        - remove_dupl_reac (bool, optional):
            Option to remove the duplicate reactions.
            True is only applied, if check_dupl_reac is also True.
            Defaults to False.
        - remove_dupl_meta (bool, optional):
            Option to remove the duplicate metabolites.
            True is only applied, if check_dupl_meta is also True.
            Defaults to False.
        - universal (str, optional):
            Path to a universal model for gapfilling.
            Defaults to None, which skips the gapfilling.
        - media_path (str, optional):
            Path to a medium config file for gapfilling.
            Defaults to None.
        - namespace (Literal['BiGG'], optional):
            Namespace to use for the model.
            Options include 'BiGG'.
            Defaults to 'BiGG'.
        - growth_threshold (float, optional):
            Growth threshold for the gapfilling.
            Defaults to 0.05.
        - iterations (int, optional):
            Number of iterations for the heuristic version of the gapfilling.
            If 0 or None is given, uses full set of reactions.
            Defaults to 3.
        - chunk_size (int, optional):
            Number of reactions to be used for gapfilling at the same time.
            If None or 0 is given, use full set, not heuristic.
            Defaults to 10000.
        - memote (bool, optional):
            Option to run memote on the cleaned model.
            Defaults to False.

    Raises:
        - ValueError: Unknown option for check_dupl_meta
        - KeyError: Missing parameter for GeneGapFiller
    """

    # -------------------
    # checking parameters
    # -------------------
    if not check_dupl_meta in ["default", "skip", "exhaustive"]:
        raise ValueError(
            "Unknown option {check_dupl_meta} for checking duplicate metabolite. Use one of: default, skip, exhaustive"
        )

    if run_gene_gapfiller:
        if not GGF_REQS.issubset(run_gene_gapfiller.keys()):
            raise KeyError(
                "At least one parameter for the GeneGapFiller is missing. Re-check your input for run_gene_gapfiller"
            )

    # -------------
    # start program
    # -------------
    logger.info(
        "\nrefinement step 2: clean-up\n################################################################################\n"
    )

    # -----------------------
    # create output directory
    # -----------------------

    try:
        Path(dir, "step2-clean-up").mkdir(parents=True, exist_ok=False)
        logger.info(f'Creating new directory {str(Path(dir,"step2-clean-up"))}')
    except FileExistsError:
        logger.info("Given directory already has required structure.")
    
    # set path for logging file
    Path(dir, "step2-clean-up", "cleanup.log").unlink(missing_ok=True)
    handler = logging.handlers.RotatingFileHandler(
        str(Path(dir, "step2-clean-up", "cleanup.log")),
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
        

    model = cobra.io.read_sbml_model(model)

    # --------------------
    # check direction
    # --------------------

    if biocyc_db:

        logger.info("\n# --------------------\n# check direction\n# --------------------")

        start = time.time()

        # check direction
        model = check_direction(model, biocyc_db)

        end = time.time()
        logger.info(f"\ttime: {end - start}s")

    # ----------
    # gapfilling
    # ----------

    logger.info("\n# ----------\n# gapfilling\n# ----------")
    start = time.time()

    # gapfilling
    ############

    # GeneGapFiller
    # -------------
    if run_gene_gapfiller:
        # load model with libsbml
        libmodel = None
        with NamedTemporaryFile(suffix=".xml", delete=False) as tmp:
            write_model_to_file(model, tmp.name)
            libmodel = load_model(tmp.name, "libsbml")
        os.remove(tmp.name)

        # run the gene gap filler
        ggf = GeneGapFiller()
        ggf.find_missing_genes(run_gene_gapfiller["gff"], libmodel)
        ggf.find_missing_reactions(
            model,
            prefix=run_gene_gapfiller["prefix"],
            type_db=run_gene_gapfiller["type_db"],
            fasta=run_gene_gapfiller["fasta"],
            dmnd_db=run_gene_gapfiller["dmnd_db"],
            map_db=run_gene_gapfiller["map_db"],
            mail=run_gene_gapfiller["mail"],
            check_NCBI=run_gene_gapfiller["check_NCBI"],
            threshold_add_reacs=run_gene_gapfiller["threshold_add_reacs"],
            outdir=Path(dir, "step2-clean-up"),
            sens=run_gene_gapfiller["sens"],
            cov=run_gene_gapfiller["cov"],
            t=run_gene_gapfiller["t"],
            pid=run_gene_gapfiller["pid"],
        )
        libmodel = ggf.fill_model(
            libmodel,
            namespace=namespace,
            idprefix=run_gene_gapfiller["prefix"],
            formula_check=run_gene_gapfiller["formula_check"],
            exclude_dna=run_gene_gapfiller["exclude_dna"],
            exclude_rna=run_gene_gapfiller["exclude_rna"],
        )

        write_model_to_file(
            libmodel, Path(dir, "step2-clean-up", "after_2ndGF_nogff.xml")
        )

        # re-load model with cobrapy
        with NamedTemporaryFile(suffix=".xml", delete=False) as tmp:
            write_model_to_file(libmodel, tmp.name)
            model = load_model(tmp.name, "cobra")
        os.remove(tmp.name)

        ggf.report(Path(dir, "step2-cleanup"))

    # gap-filling via COBRApy medium
    # ------------------------------

    # construct a list of media
    media_list = []

    # load media from config file
    if media_path:
        media_list = load_media(media_path)

    #   separate option for cobra gapfilling
    if len(media_list) > 0:
        # load universal model
        universal_model = load_model(universal, "cobra")
        # run gapfilling
        model = multiple_cobra_gapfill(
            model,
            universal_model,
            media_list,
            namespace,
            iterations=iterations,
            chunk_size=chunk_size,
            growth_threshold=growth_threshold,
        )

    end = time.time()
    logger.info(f"\ttime: {end - start}s")

    # -----------------
    # resove duplicates
    # -----------------
    logger.info("\n# -----------------\n# resolve duplicates\n# -----------------")
    start = time.time()

    model = resolve_duplicates(
        model,
        check_reac=check_dupl_reac,
        check_meta=check_dupl_meta,
        remove_unused_meta=remove_unused_meta,
        remove_dupl_reac=remove_dupl_reac,
        replace_dupl_meta=remove_dupl_meta,
    )

    end = time.time()
    logger.info(f"\ttime: {end - start}s")

    # ---------------------
    # dead ends and orphans
    # ---------------------
    # (currently no removal of dead ends and orphans as they may be interesting
    # for manual curation)

    # ----------
    # save model
    # ----------
    name = f"{model.id}_clean"
    cobra.io.write_sbml_model(model, Path(dir, "step2-clean-up", name + ".xml"))

    # ---------------------------------
    # assess model quality using memote
    # ---------------------------------

    if memote:
        memote_path = str(Path(dir, "step2-clean-up", name + ".html"))
        run_memote(model, "html", return_res=False, save_res=memote_path, verbose=True)


# Part 3: annotation
# ------------------


def annotate(
    model: str,
    dir: str,
    kegg_viaEC: bool = False,
    kegg_viaRC: bool = False,
    memote: bool = False,
):
    """Further annotate a given model.

    Currently can add annotations for:

    - SBO using SBOannotator
    - KEGG.reaction -> KEGG.pathway

    Args:
        - model (str):
            Path to the model.
        - dir (str):
            Path to the output directory.
        - kegg_viaEC (bool, optional):
            Option to search for KEGG pathway ID using the
            EC number if previous searches were unsuccesful.
            Defaults to False.
        - kegg_viaRC (bool, optional):
            Option to search for KEGG pathway ID using the
            reaction class if previous searches were unsuccesful.
            Defaults to False.
        - memote (bool, optional):
            Optionally run memote after the annotation.
            Defaults to False.
    """

    logger.info(
        "\nrefinement step 3: annotation\n################################################################################\n"
    )

    # -----------------------
    # create output directory
    # -----------------------

    try:
        Path(dir, "step3-annotation").mkdir(parents=True, exist_ok=False)
        logger.info(f'Creating new directory {str(Path(dir,"step3-annotation"))}')
    except FileExistsError:
        logger.info("Given directory already has required structure.")
        
    # set path for logging file
    Path(dir, "step3-annotation", "annotation.log").unlink(missing_ok=True)
    handler = logging.handlers.RotatingFileHandler(
        str(Path(dir, "step3-annotation", "annotation.log")),
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

    # --------------
    # Load the model
    # --------------

    model = load_model(model, "libsbml")

    # ------------------
    # Polish annotations
    # ------------------

    logger.info(
        "\n# ----------------------------------\n# polish annotations\n# ----------------------------------"
    )
    start = time.time()

    model = polish_annotations(
        model,
        True,
        str(Path(dir, "step3-annotation")),
    )

    end = time.time()
    logger.info(f"\ttime: {end - start}s")

    # ------------------
    # add SBO annotation
    # ------------------

    logger.info("\n# ------------------\n# add SBO annotation\n# ------------------")

    start = time.time()

    model = run_SBOannotator(model)
    write_model_to_file(
        model, str(Path(dir, "step3-annotation", model.getId() + "_SBOannotated.xml"))
    )

    end = time.time()
    logger.info(f"\ttime: {end - start}s")

    # ----------------
    # add KEGG pathway
    # ----------------
    logger.info("\n# ----------------\n# add KEGG pathway\n# ----------------")

    start = time.time()

    libmodel, nokegg = set_kegg_pathways(
        str(Path(dir, "step3-annotation", model.getId() + "_SBOannotated.xml")),
        viaEC=kegg_viaEC,
        viaRC=kegg_viaRC,
    )
    write_model_to_file(
        libmodel,
        str(Path(dir, "step3-annotation", model.getId() + "_keggpathways.xml")),
    )
    with open(str(Path(dir, "step3-annotation", "reac_no_keggORec")), "w") as f:
        for line in nokegg:
            f.write(f"{line}\n")

    end = time.time()
    logger.info(f"\ttime: {end - start}s")

    # ---------------------------------
    # assess model quality using memote
    # ---------------------------------

    if memote:
        start = time.time()
        model = load_model(
            str(
                Path(dir, "step3-annotation", model.getId() + "_keggpathways.xml"),
                "cobra",
            )
        )
        name = model.id
        memote_path = str(Path(dir, "step3-annotation", name + "_annotated.html"))
        run_memote(model, "html", return_res=False, save_res=memote_path, verbose=True)
        end = time.time()
        logger.info(f"\ttotal time: {end - start}s")


# Part 4: smoothing
# -----------------


def smooth(
    genome: str,
    model: str,
    dir: str,
    mcc="skip",
    egc_solver: None | Literal["greedy"] = None,
    namespace: Literal["BiGG"] = "BiGG",
    dna_weight_frac=0.023,
    ion_weight_frac=0.05,
    memote=False,
):
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

    logger.info(
        "\nrefinement step 4: smoothing\n################################################################################\n"
    )

    # -----------------------
    # create output directory
    # -----------------------

    try:
        Path(dir, "step4-smoothing").mkdir(parents=True, exist_ok=False)
        logger.info(f'Creating new directory {str(Path(dir,"step4-smoothing"))}')
    except FileExistsError:
        logger.info("Given directory already has required structure.")

    try:
        Path(dir, "manual_curation").mkdir(parents=True, exist_ok=False)
        logger.info(f'Creating new directory {str(Path(dir,"manual_curation"))}')
    except FileExistsError:
        logger.info("Given directory already has required structure.")
        
    # set path for logging file
    Path(dir, "step4-smoothing", "smoothing.log").unlink(missing_ok=True)
    handler = logging.handlers.RotatingFileHandler(
        str(Path(dir, "step4-smoothing", "smoothing.log")),
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

    # ---------
    # load data
    # ---------
    model = cobra.io.read_sbml_model(model)

    # ---------------
    # mass and charge
    # ---------------

    if mcc == "apply":
        logger.info(
            "\n# ----------------------------------\n# mass and charge curation (applied)\n# ----------------------------------"
        )
        start = time.time()
        model = perform_mcc(model, Path(dir, "step4-smoothing"))
        end = time.time()
        logger.info(f"\ttime: {end - start}s")
    elif mcc == "extra":
        logger.info(
            "\n# --------------------------------\n# mass and charge curation (extra)\n# --------------------------------"
        )
        start = time.time()
        model = perform_mcc(model, Path(dir, "manual_curation"), False)
        end = time.time()
        logger.info(f"\ttime: {end - start}s")
    elif mcc == "skip":
        logger.info(
            "\n# ------------------------\n# mass and charge curation\n# ------------------------\n\tskipped"
        )
    else:
        warnings.warn(
            f"Unknown option {mcc} for Mass and Charge Curation. Usage of MCC will be skipped."
        )

    # ----------------------------------
    # check for energy generating cycles
    # ----------------------------------
    logger.info(
        "\n# ---------------------------------------------\n# # check for energy generating cycles\n# ---------------------------------------------"
    )
    start = time.time()

    match egc_solver:
        # greedy solver
        case "greedy":
            logger.info("Using GreedyEGCSolver...")
            solver = egcs.GreedyEGCSolver()
            results = solver.solve_egcs(
                model, namespace=namespace
            )  # automatically uses c,e as compartments
            if results:
                for k, v in results.items():
                    logger.info(f"\t{k}: {v}")

        # no solver = EGCs will only be reported
        case _:
            solver = egcs.EGCSolver()
            logger.info(f"\tFound EGCs:\n")
            logger.info(
                f"\t{solver.find_egcs(model,with_reacs=True,namespace=namespace)}"
            )  # automatically uses c,e as compartments

    end = time.time()
    logger.info(f"\ttime: {end - start}s")

    # --------------------------
    # biomass objective function
    # --------------------------
    # adjust the BOF to the current genome

    logger.info("\n# ----------\n# adjust BOF\n# ----------")
    start = time.time()

    with tempfile.NamedTemporaryFile(suffix=".xml", delete=False) as temp_model:
        # generate an up-to-date model xml-file
        cobra.io.write_sbml_model(model, temp_model.name)
        # update BOF
        pos_bofs = test_biomass_presence(model)
        if pos_bofs:
            model.reactions.get_by_id(pos_bofs[0]).reaction = adjust_BOF(
                genome, temp_model.name, model, dna_weight_frac, ion_weight_frac
            )
            # optimise BOF(s)
            model = check_normalise_biomass(model)
        else:
            # create new BOF
            bof_reac = Reaction("Biomass_BOFdat")
            bof_reac.name = "Biomass objective function created by BOFdat"
            model.add_reactions([bof_reac])
            model.reactions.get_by_id(bof_reac).reaction = adjust_BOF(
                genome, temp_model.name, model, dna_weight_frac, ion_weight_frac
            )

            # optimise BOF(s)
            model = check_normalise_biomass(model)
    os.remove(temp_model.name)

    end = time.time()
    logger.info(f"\ttime: {end - start}s")

    # ----------------
    # save final model
    # ----------------
    logger.info("\n# ----------\n# save model\n# ----------")
    model_name = f"{model.id}_smooth"
    outname = Path(dir, "step4-smoothing", model_name + ".xml")
    logger.info(f"\tsaving to: {outname}")
    cobra.io.write_sbml_model(model, outname)

    # ---------------------------------
    # assess model quality using memote
    # ---------------------------------

    if memote:
        memote_path = str(Path(dir, "step4-smoothing", model_name + ".html"))
        run_memote(model, "html", return_res=False, save_res=memote_path, verbose=True)
