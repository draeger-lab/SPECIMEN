"""Generate a draft model from a template.
"""
__author__ = 'Carolin Brune'
################################################################################
# requirements
################################################################################

from pathlib import Path
import os.path
import time
import pandas as pd
import cobra
import numpy as np
import subprocess
import warnings

from refinegems.medium import load_medium_from_db, medium_to_model
from refinegems.io import load_model_cobra

from ..util import cobra_models

# further required programs:
#        - MEMOTE, tested with version 0.13.0

################################################################################
# functions
################################################################################

def pid_filter(data: pd.DataFrame, pid: float):
    """Filter the data based on PID threshold. Entries above the given value are kept.

    :param data: The data from bidirectional_blast.py containing at least a 'PID' column.
    :type  data: pd.DataFrame
    :param pid:  PID threshold value.
    :type  pid:  double, should be given in percentage e.g. 80.0.
    :returns:    The filtered data.
    :rtype:      pd.DataFrame
    """

    data['homolog_found'] = np.where(data['PID'] > pid, 1, 0)
    counts = data['homolog_found'].value_counts()
    print(F'\t{counts[1]} set as 1 (= homologs), {counts[0]} set as 0 ')

    return data

# @TODO
def edit_template_identifiers(data, edit):
    """Edit the subject IDs to fit the gene IDs of the template model.
    Requires further extention, if needed edits are not included.

    :param data: the data frame containing the bidirectional blastp best hits information.
    :type  data: pd.DataFrame
    :param edit: type of edit that is to be performed
    :type  edit: string, currently implemented: no, dot-to-underscore
    :returns:    the edited data
    :rtype:      pd.DataFrame
    """

    match edit:
        case 'no':
            return data
        case 'dot-to-underscore':
            data['subject_ID'] = [x.replace('.','_') for x in data['subject_ID']]
            return data
        # ...........
        # @TODO
        # extendable
        # ...........
    return data


def remove_absent_genes(model, genes):
    """Remove a list of genes from a given model.
    Note: genes that are not found in the model are skipped.

    :param model: A model to delete genes from. A copy will be created before deleting.
    :type  model: cobra.Model
    :param genes: Gene identifiers of genes that should be deleted.
    :type  genes: list, of strings
    :returns:     A new model with the given genes deleted, if found in the original model.
    :rtype:       cobra.Model
    """

    print('\tremove absent (low PID) genes')
    print('\t...................................................................')

    genes_to_delete = 0
    essential = 0
    not_found = 0
    remove = False

    modelCopy = model.copy()
    for g in genes:
        try:
            test = modelCopy.genes.get_by_id(g)
            # check, if gene is essential
            with modelCopy as variant:
                test.knock_out()
                variant.optimize()
                # set gene for deletion if model still works without it
                if variant.objective.value > cobra_models.MIN_GROWTH_RATE:
                    remove = True
                    genes_to_delete += 1
                # keep gene, if model stops growing
                else:
                    essential += 1
                    continue
            # remove gene
            if remove:
                cobra.manipulation.delete.remove_genes(modelCopy, [g], remove_reactions=True)
                remove = False
        except KeyError:
            not_found += 1

    print(F'\tnumber of deleted genes: {genes_to_delete}')
    print(F'\tnumber of essential genes (not deleted): {essential}')
    print(F'\tnumber of genes not found in the model: {not_found}')
    print('\t...................................................................')

    return modelCopy


def rename_found_homologs(draft, bbh):
    """Rename the genes in the model correnspondingly to the homologous ones found in the query.

    :param draft: The draft model with the to-be-renamed genes.
    :type  draft: cobra.Model
    :param bbh:   The table from the bidirectional_blast.py script containing the bidirectional blastp best hits information
    :type  bbh:   pd.DataFrame
    :returns:     The draft model with renamed genes.
    :rtype:       cobra.Model
    """

    print('\trename found homologs')
    print('\t...................................................................')
    # extract found homologs
    present_s = bbh[bbh['homolog_found'] == 1]['subject_ID'].tolist()
    present_q = bbh[bbh['homolog_found'] == 1]['query_ID'].tolist()

    print(F'\ttotal number of found homologs: {len(present_q)}')

    # map the new names to the model, if the original gene is found
    name_mapping = dict(zip(present_s, present_q))
    cobra.manipulation.modify.rename_genes(draft, name_mapping)
    renamed = len([x for x in present_q if x in draft.genes])
    print(F'\tnumber of homologs found and renamed in model: {renamed}')

    # identify skipped genes of the query
    skipped_genes = [_ for _ in present_q if _ not in name_mapping.values()]
    remap_skipped = {}
    for skp in skipped_genes:
        key = name_mapping[bbh.loc[bbh['query_ID'] == skp, 'subject_ID'].values[0]]
        if key in remap_skipped.keys():
            remap_skipped[key].append(skp)
        else:
            remap_skipped[key] = [skp]

    # create new genes entry from the old ones for additional homologous genes
    skp_counter = 0
    for k,v in remap_skipped.items():
        if k in [_.id for _ in draft.genes]:
            for ident in v:
                skp_counter += 1
                # add the gene as an alternative in the gene reaction rules of the model
                for r in draft.genes.get_by_id(k).reactions:
                    r.gene_reaction_rule += F' or {ident}'
                # copy the annotations of the homologous one for better model quality
                draft.genes.get_by_id(ident).annotation = draft.genes.get_by_id(k).annotation

    print(F'\tnumber of additional homologs found and added to model: {skp_counter}')

    print('\t...................................................................')

    return draft


def check_unchanged(draft, bbh):
    """Check the genes names (more correctly, the IDs) for still existing original col_names.
    Depending on the case, decide if to keep or remove them.

    :param draft:            The draft model currently in the making.
    :type  draft:            cobra.Model
    :param bbh:              The table from the bidirectional_blast.py script containing the bidirectional blastp best hits information
    :type  bbh:              pd.DataFrame
    :returns:                The edited draft model.
    :rtype:                  cobra.Model
    """

    print('\tcheck not renamed genes')
    print('\t...................................................................')

    query_ids = bbh['query_ID'].tolist()
    not_renamed = [x for x in draft.genes if not x.id in query_ids]
    print(F'\tnumber of not renamed genes: {len(not_renamed)}')

    to_delete = 0
    essential_counter = 0
    remove = False

    for g in not_renamed:
        with draft as model:
            g.knock_out()
            model.optimize()
            if model.slim_optimize() > cobra_models.MIN_GROWTH_RATE:
                # set for deletion
                remove = True
                to_delete += 1
            else:
                # keep
                essential_counter += 1

        # delete gene
        if remove:
            cobra.manipulation.delete.remove_genes(draft, [g], remove_reactions=True)
            remove = False

    print(F'\tremoved {to_delete} non-essential genes')
    print(F'\tkept {essential_counter} essential genes')
    print('\t...................................................................')

    return draft


def gen_draft_model(model, bbh, name, dir, edit, medium='default', namespace='BiGG'):
    """Generate a draft model from a template model and the results of a bidirectional blastp (blast best hits) table
    and save it as a new model.

    :param model:            The template model.
    :type  model:            cobra.Model
    :param bbh:              The bidirectional blastp best hits table.
    :type  bbh:              pd.DataFrame
    :param name:             The name of the generated model.
    :type  name:             string
    :param dir:              Path to the directory to save the new model in.
    :type  dir:              string
    :param edit:             Allows for editing the subject_ID od the bbh table to be edited
    :type  edit:             string, currently implemented are options no and dot-to-underscore
    :param medium: Name of the to be loaded from a database or 'default' = the one
        from the template model.
    :type medium: string
    :param namespace:        Namespace of the model.
    :type  namespace:        Literal['BiGG'], optional
    """

    match medium:
        # use the medium from the template
        case 'default':
            pass
        # set all exchanges open + as medium
        case 'exchanges':
            # @NOTE: gapfilling using cobra becomes not feasible using this option
            model.medium = {_.id:1000.0 for _ in model.exchanges}
        # use a medium from the refinegems database
        case str():
            new_m = load_medium_from_db(medium)
            medium_to_model(model, new_m, namespace, double_o2=False, add=True)
        case _:
            warnings.warn('Unknown or incomplete input for setting the medium. Using the one from the template.')

    # delete absent genes, including associated reactions
    bbh = edit_template_identifiers(bbh, edit)
    absent = set(bbh[bbh['homolog_found'] == 0]['subject_ID'].tolist())

    draft = remove_absent_genes(model, absent)

    # rename genes as per the naming of the new model (if possible)
    draft = rename_found_homologs(draft, bbh)

    # check genes, that have not gotten a new ID assigned
    draft = check_unchanged(draft, bbh)

    # rename compartments to the standart
    draft = cobra_models.resolve_compartment_names(draft)

    # for each object, save a note that it was added during draft construction
    for r in draft.reactions:
        r.notes['creation'] = 'via template'
    for g in draft.genes:
        g.notes['creation'] = 'via template'
    for m in draft.metabolites:
        m.notes['creation'] = 'via template'

    # save draft model
    draft.id = name
    cobra.io.write_sbml_model(draft, F'{dir}{name}_draft.xml')


def run(template, bpbbh, dir, edit_names='no', pid=80.0, name=None, medium='default', namespace='BiGG', memote=False):
    """Generate a draft model from a blastp best hits tsv file and a template model.

    :param template: Path to the file containing the template model.
    :type template: string
    :param bpbbh: Path to the blastp bidirectional best hits.
    :type bpbbh: string
    :param dir: Path to output directory.
    :type dir: string
    :param edit_names: Edit the identifier of the FASTA files to match the names in the model.
        Default is no. Currently implemented edits: dot-to-underscore'.
    :type edit_names: string
    :param pid: Threshold value for determining, if a gene is counted as present or absent.
        The default is 80.0 (percentage).
    :type pid: float
    :param name: Name of the output.
    :type name: string, optional
    :param medium: Set the medium for the new model. If not set, will use the one from the template.
        If given the keyword "exchanges", will open all exchange reaction and use them as the medium.
        If given together with db_path, the corresponding medium is loaded.
    :type medium: string
    :param namespace: Namespace of the model.
    :type  namespace: Literal['BiGG'], optional
    :param memote: Option to run memote after creating the draft model.
        Default is False (not running memote).
    :type memote: bool

    :raises: :class:`ValueError`: 'Edit_names value {edit_names} not in list of allowed values: no, dot-to-underscore'
    :raises: :class:`ValueError`: 'Model type {model_type} not in list of allowed values: sbml, json, matlab'
    """

    total_time_s = time.time()

    # -----------
    # check input
    # -----------

    if edit_names not in ['no','dot-to-underscore']:
        raise ValueError(F'Edit_names value {edit_names} not in list of allowed values: no, dot-to-underscore')

    if name is None:
        name = "_".join(os.path.splitext(bpbbh)[0].split('_', 2)[:2])

    # -------------
    # start program
    # -------------

    print('\ngenerate draft model\n################################################################################\n')

    try:
        Path(dir).mkdir(parents=True, exist_ok=False)
        print(F'Creating new directory {dir}')
    except FileExistsError:
        print('Given directory already exists.')

    bbh_data = pd.read_csv(bpbbh, sep='\t')
    template_model = load_model_cobra(template)
    # if possible, set growth function as model objective
    growth_objfunc = cobra_models.find_growth_obj_func(template_model)
    template_model.objective = growth_objfunc

    # -----------------------------
    # determine presence or absence
    # -----------------------------

    print('\n# ------------------\n# filter by PID\n# ------------------')
    start = time.time()
    pid_filter(bbh_data, pid)
    end = time.time()
    print(F'\ttotal time: {end - start}s')

    # --------------------
    # generate draft model
    # --------------------

    print('\n# --------------------\n# generate draft model\n# --------------------')
    start = time.time()
    gen_draft_model(template_model, bbh_data, name, dir, edit_names, medium=medium, namespace=namespace)
    end = time.time()
    print(F'\ttotal time: {end - start}s')

    # -------------------
    # analyse with MEMOTE
    # -------------------
    if memote:
        print('\n# -------------------\n# analyse with MEMOTE\n# -------------------')
        start = time.time()
        draft_path = F'{dir}{name}_draft.xml'.replace(" ", "\ ")
        memote_path = F'{dir}{name}_draft.html'.replace(" ", "\ ")
        subprocess.run([F'memote report snapshot --filename {memote_path} {draft_path}'], shell=True)
        end = time.time()
        print(F'\ttotal time: {end - start}s')

    total_time_e = time.time()
    print(F'total runtime: {total_time_e-total_time_s}')
