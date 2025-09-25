"""Generate a draft model from a template model.

The basic idea has been adapted from Norsigian et al. (2020).
"""

__author__ = "Carolin Brune"
################################################################################
# requirements
################################################################################

import cobra
import importlib.metadata
import logging
import numpy as np
import os.path
import pandas as pd
import time

from pathlib import Path
from typing import Literal, Union

from refinegems.classes.medium import load_medium_from_db, medium_to_model
from refinegems.utility.io import load_model
from refinegems.utility.entities import (
    resolve_compartment_names,
    remove_non_essential_genes,
)
from refinegems.utility.util import test_biomass_presence
from refinegems.utility.connections import run_memote
from refinegems.curation.curate import fix_reac_bounds
from refinegems.utility.cvterms import add_cv_term_reactions

from refinegems.utility.util import MIN_GROWTH_THRESHOLD

# further required programs:
#        - MEMOTE, tested with version 0.13.0+

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


def pid_filter(data: pd.DataFrame, pid: float) -> pd.DataFrame:
    """Filter the data based on PID threshold. Entries above the given value are 
    retained.

    Args:
        - data (pd.DataFrame):
            The data from teh previous step (see :py:mod:`~specimen.hqtb.core.bidirectional_blast`) 
            containing at least a 'PID' column.
        - pid (float):
            PID threshold value, given in percentage e.g. 80.0.

    Returns:
        pd.DataFrame:
            The filtered data.
    """
    # get logger
    logger = logging.getLogger(__name__ + "-intern")

    # filter data
    data["homolog_found"] = np.where(data["PID"] > pid, 1, 0)
    counts = data["homolog_found"].value_counts()
    logger.info(f"\t{counts[1]} set as 1 (= homologs), {counts[0]} set as 0 ")

    return data


def edit_template_identifiers(
    data: pd.DataFrame, edit: Literal["no", "dot-to-underscore"]
) -> pd.DataFrame:
    """Edit the subject IDs to fit the gene IDs of the template model.
    Requires further extention, if needed edits are not included.

    Args:
        - data (pd.DataFrame):
            The data frame containing the bidirectional blastp best hits information.
        - edit (Literal['no','dot-to-underscore']):
            Type of edit to perform.
            Currently possible options: no, dot-to-underscore.

    Returns:
        pd.DataFrame:
            The (un)edited DataFrame.
    """
    # get logger
    logger = logging.getLogger(__name__ + "-intern")

    match edit:
        case "no":
            return data
        case "dot-to-underscore":
            data["subject_ID"] = [x.replace(".", "_") for x in data["subject_ID"]]
            return data
        case _:
            mes = "Unknown option for parameter edit. Nothing will be edited. If you need another option, contact the developers."
            logger.warning(mes)
            return data


def remove_absent_genes(model: cobra.Model, genes: list[str]) -> cobra.Model:
    """Remove a list of genes from a given model.
    
    .. note:: 
        Genes that are not found in the model are skipped.

    Args:
        - model (cobra.Model):
            A template model to delete genes from.
            A copy will be created before deleting.
        - genes (list[str]):
            Gene identifiers of genes that should be deleted.

    Returns:
        cobra.Model:
            A new model with the given genes deleted, if found in the original model.
    """
    # get logger
    logger = logging.getLogger(__name__ + "-intern")

    logger.info("\tremove absent (low PID) genes")
    logger.info("\t...................................................................")

    genes_to_delete = 0
    essential = 0
    incomplex = 0
    not_found = 0
    remove = False

    # @NOTE
    modelCopy = model.copy() # currently not compatible with 3.13
    # modelCopy = copy.deepcopy(model) # also not working with 3.13 
    # there is a - supposed - fix on cobrapy, but not yet released
    for g in genes:
        try:
            test = modelCopy.genes.get_by_id(g)
            # check, if gene is essential
            with modelCopy as variant:
                test.knock_out()
                variant.optimize()
                # set gene for deletion if model still works without it
                if variant.objective.value > MIN_GROWTH_THRESHOLD:
                    complexctr = False
                    # check if gene is part of a partially remapped complex 
                    # if yes, keep it regardless of essentiality
                    for r in modelCopy.reactions:
                        reac_genes = [_.id for _ in list(r.genes)]
                        if reac_genes != 0 and ('and '+g in r.gene_reaction_rule or g+' and' in r.gene_reaction_rule) and not set(reac_genes).issubset(genes):
                            logger.info(f'Keeping gene {g}, as it is part of a complex, where not all genes should be deleted.')
                            complexctr = True 
                            break
                    if complexctr:
                        incomplex += 1
                    else:
                        remove = True
                        genes_to_delete += 1
                    
                # keep gene, if model stops growing
                else:
                    essential += 1
                    continue
                
            # remove gene
            if remove:
                cobra.manipulation.delete.remove_genes(
                    modelCopy, [g], remove_reactions=True
                )
                remove = False
        except KeyError as e:
            logger.warning(f"Gene {g} could not be found.")
            not_found += 1

    logger.info(f"\tnumber of deleted genes: {genes_to_delete}")
    logger.info(f"\tnumber of essential genes (not deleted): {essential}")
    logger.info(f"\tnumber of genes not found in the model: {not_found}")
    logger.info(f"\tnumber of genes not deleted due to them being part of a complex: {incomplex}")
    logger.info("\t...................................................................")

    return modelCopy


def rename_found_homologs(draft: cobra.Model, bbh: pd.DataFrame) -> cobra.Model:
    """Rename the genes in the model correnspondingly to the homologous ones found in the query.

    Args:
        - draft (cobra.Model):
            The draft model with the to-be-renamed genes.
        - bbh (pd.DataFrame):
            The table from :py:func:`~specimen.hqtb.core.bidirectional_blast.run` containing 
            the bidirectional blastp best hits information

    Returns:
        cobra.Model:
            The draft model with renamed genes.
    """
    # get logger
    logger = logging.getLogger(__name__ + "-intern")

    logger.info("\trename found homologs")
    logger.info("\t...................................................................")
    # extract found homologs
    present_s = bbh[bbh["homolog_found"] == 1]["subject_ID"].tolist()
    present_q = bbh[bbh["homolog_found"] == 1]["query_ID"].tolist()

    logger.info(f"\ttotal number of found homologs: {len(present_q)}")

    # map the new names to the model, if the original gene is found
    name_mapping = dict(zip(present_s, present_q))
    cobra.manipulation.modify.rename_genes(draft, name_mapping)
    renamed = len([x for x in present_q if x in draft.genes])
    logger.info(f"\tnumber of homologs found and renamed in model: {renamed}")

    # identify skipped genes of the query
    skipped_genes = [_ for _ in present_q if _ not in name_mapping.values()]
    remap_skipped = {}
    
    for skp in skipped_genes:
        key = name_mapping[bbh.loc[bbh["query_ID"] == skp, "subject_ID"].values[0]]
        if key in remap_skipped.keys():
            remap_skipped[key].append(skp)
        else:
            remap_skipped[key] = [skp]

    # create new genes entry from the old ones for additional homologous genes
    skp_counter = 0
    for k, v in remap_skipped.items():
        if k in [_.id for _ in draft.genes]:
            for ident in v:
                skp_counter += 1
                # add the gene as an alternative in the gene reaction rules of the model
                for r in draft.genes.get_by_id(k).reactions:
                    rgpr = r.gene_reaction_rule
                    if k+' and' in rgpr or 'and '+k in  rgpr:
                        r.gene_reaction_rule = rgpr.replace(k,f'({k} or {ident})') 
                    else:
                        r.gene_reaction_rule = rgpr.replace(k,f'{k} or {ident}') 
                    # copy the annotations of the homologous one for better model quality
                    draft.genes.get_by_id(ident).annotation = draft.genes.get_by_id(k).annotation

    logger.info(
        f"\tnumber of additional homologs found and added to model: {skp_counter}"
    )

    logger.info("\t...................................................................")

    return draft


def check_unchanged(draft: cobra.Model, bbh: pd.DataFrame) -> cobra.Model:
    """Check the genes names (more correctly, the IDs) for still existing original col_names.
    Depending on the case, decide if to keep or remove them.

    Args:
        - draft (cobra.Model):
            The draft model currently in the making.
        - bbh (pd.DataFrame):
            The table from :py:func:`~specimen.hqtb.core.bidirectional_blast.run` 
            containing the bidirectional blastp best hits information.

    Returns:
        cobra.Model:
            The model after the check and possible removal of genes.
    """
    # get logger
    logger = logging.getLogger(__name__ + "-intern")

    logger.info("\tcheck not renamed genes")
    logger.info("\t...................................................................")

    query_ids = bbh["query_ID"].tolist()
    not_renamed = [x for x in draft.genes if not x.id in query_ids]
    logger.info(f"\tnumber of not renamed genes: {len(not_renamed)}")

    to_delete, essential_counter, in_complex = remove_non_essential_genes(
        draft, genes_to_check=not_renamed
    )

    logger.info(f"\tremoved {to_delete} non-essential genes")
    logger.info(f"\tkept {essential_counter} essential genes")
    logger.info(f"\tkept {in_complex} (non-essential) genes found in complexes")
    logger.info("\t...................................................................")

    return draft


def gen_draft_model(
    model: cobra.Model,
    bbh: pd.DataFrame,
    name: str,
    dir: str,
    edit: Literal["no", "dot-to-underscore"],
    medium: str = "default",
    namespace: Literal["BiGG"] = "BiGG",
) -> cobra.Model:
    """Generate a draft model from a template model and the results of a bidirectional blastp (blast best hits) table
    and save it as a new model.

    Args:
        - model (cobra.Model):
            The template model.
        - bbh (pd.DataFrame):
            The bidirectional blastp best hits table.
        - name (str):
            Name of the newly generated model.
        - dir (str):
            Path to the directory to save the new model in.
        - edit (Literal['no','dot-to-underscore'):
            Type of edit to perform.
            Currently possible options: no, dot-to-underscore.
        - medium (str, optional):
            Name of the to be loaded from the refineGEMs database or 'default' = the one
            from the template model. If given the keyword 'exchanges', will use all exchange reactions in the model as a medium.
            Defaults to 'default'.
        - namespace (Literal['BiGG'], optional):
            Namespace of the model.
            Defaults to 'BiGG'.

    Returns:
        cobra.Model:
            The generated draft model.
    """
    # get logger
    logger = logging.getLogger(__name__ + "-intern")

    match medium:
        # use the medium from the template
        case "default":
            pass
        # set all exchanges open + as medium
        case "exchanges":
            logger.warning(
                "Medium set to exchanges. If gap filling with COBRApy is enabled, there is a high possibility, the program will crash with an infeasible error."
            )
            model.medium = {_.id: 1000.0 for _ in model.exchanges}
        # use a medium from the refinegems database
        case str():
            new_m = load_medium_from_db(medium)
            medium_to_model(model, new_m, namespace, double_o2=False, add=True)
        case _:
            logger.warning(
                "Unknown or incomplete input for setting the medium. Using the one from the template."
            )

    # delete absent genes, including associated reactions
    bbh = edit_template_identifiers(bbh, edit)
    absent = set(bbh[bbh["homolog_found"] == 0]["subject_ID"].tolist())

    draft = remove_absent_genes(model, absent)

    # rename genes as per the naming of the new model (if possible)
    draft = rename_found_homologs(draft, bbh)

    # check genes, that have not gotten a new ID assigned
    draft = check_unchanged(draft, bbh)

    # rename compartments to the standard
    resolve_compartment_names(draft)
    # smooth out potential issues with reaction bounds
    fix_reac_bounds(draft)

    # for each object, save a note that it was added during draft construction
    for r in draft.reactions:
        r.notes["creation"] = "via template"
        add_cv_term_reactions("0007482", "ECO", r)
    for g in draft.genes:
        g.notes["creation"] = "via template"
    for m in draft.metabolites:
        m.notes["creation"] = "via template"

    # save draft model
    old_desc = draft.id if draft.id else (draft.name if draft.name else "Unknown template")
    draft.id = name
    draft.notes["Template model"] = old_desc
    draft.notes["Description"] = (
        f'This model was created with SPECIMEN version {importlib.metadata.version("specimen")}'
    )
    draft.name = f"Genome-scale metabolic model {name}"
    cobra.io.write_sbml_model(draft, Path(dir, name + "_draft.xml"))

    return draft


def run(
    template: str,
    bpbbh: str,
    dir: str,
    edit_names: Literal["no", "dot-to-underscore"] = "no",
    pid: float = 80.0,
    name: Union[str, None] = None,
    medium: str = "default",
    namespace: str = "BiGG",
    memote: bool = False,
):
    """Generate a draft model from a blastp best hits tsv file and a template model.

    Args:
        - template (str):
            Path to the file containing the template model.
        - bpbbh (str):
            Path to the blastp bidirectional best hits.
        - dir (str):
            Path to output directory.
        - edit_names (Literal['no','dot-to-underscore', optional):
            Type of edit to perform.
            Currently possible options: no, dot-to-underscore.
            Defaults to 'no'.
        - pid (float, optional):
            Threshold value for determining, if a gene is counted as present or absent.
            Given in percentage, e.g. 80.0 = 80%.
            Defaults to 80.0.
        - name (Union[str,None], optional):
            Name of the output model.
            If not given, takes name from filename.
            Defaults to None.
        - medium (str, optional):
            Name of the medium to be loaded from the refineGEMs database or 'default' = the one
            from the template model. If given the keyword 'exchanges', will use all exchange reactions in the model as a medium.
            Defaults to 'default'.
        - namespace (str, optional):
            Namespace of the model.
            Defaults to 'BiGG'.
        - memote (bool, optional):
            Option to run memote after creating the draft model.
            Defaults to False.

    Raises:
        - ValueError: 'Edit_names value not in list of allowed values: no, dot-to-underscore'
    """
    total_time_s = time.time()

    # -----------------------
    # create output directory
    # -----------------------

    try:
        Path(dir).mkdir(parents=True, exist_ok=False)
        genlogger.info(f"Creating new directory {dir}")
    except FileExistsError:
        genlogger.info(f"Given directory already exists: {dir}")

    # -----------------
    # fine tune logging
    # -----------------
    # interal logging
    Path(dir, "generate_draft_model.log").unlink(missing_ok=True)
    handler = logging.handlers.RotatingFileHandler(
        str(Path(dir, "generate_draft_model.log")),
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

    # -----------
    # check input
    # -----------
    if edit_names not in ["no", "dot-to-underscore"]:
        raise ValueError(
            f"Edit_names value {edit_names} not in list of allowed values: no, dot-to-underscore"
        )

    if name is None:
        name = "_".join(os.path.splitext(bpbbh)[0].split("_", 2)[:2])

    # -------------
    # start program
    # -------------

    logger.info(
        "\ngenerate draft model\n################################################################################\n"
    )

    bbh_data = pd.read_csv(bpbbh, sep="\t")
    # ensure, IDs are seen as strings
    bbh_data["query_ID"] = bbh_data["query_ID"].astype(str)
    bbh_data["subject_ID"] = bbh_data["subject_ID"].astype(str)
    # load template
    template_model = load_model(template, "cobra")
    # if possible, set growth function as model objective
    growth_objfunc = test_biomass_presence(template_model)
    if len(growth_objfunc) == 1:
        template_model.objective = growth_objfunc[0]
    elif len(growth_objfunc) > 1:
        mes = f"Multiple BOF detected. Choosing the following: {growth_objfunc[0]}"
        logger.warning(mes)
        template_model.objective = growth_objfunc[0]
    else:
        mes = f"No BOF detected. Can lead to problems downstream the SPECIMEN pipeline."
        logger.warning(mes)

    # -----------------------------
    # determine presence or absence
    # -----------------------------

    logger.info("\n# ------------------\n# filter by PID\n# ------------------")
    start = time.time()
    pid_filter(bbh_data, pid)
    end = time.time()
    logger.info(f"\ttotal time: {end - start}s")

    # --------------------
    # generate draft model
    # --------------------

    logger.info(
        "\n# --------------------\n# generate draft model\n# --------------------"
    )
    start = time.time()
    draft = gen_draft_model(
        template_model,
        bbh_data,
        name,
        dir,
        edit_names,
        medium=medium,
        namespace=namespace,
    )
    end = time.time()
    logger.info(f"\ttotal time: {end - start}s")

    # -------------------
    # analyse with MEMOTE
    # -------------------

    if memote:
        memote_path = str(Path(dir, name + ".html"))
        run_memote(draft, "html", return_res=False, save_res=memote_path, verbose=True)

    total_time_e = time.time()
    logger.info(f"total runtime: {total_time_e-total_time_s}")

    # restore logging behaviour
    cobralogger.handlers.clear()
    cobralogger.propagate = False    
