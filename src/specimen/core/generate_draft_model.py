"""Generate a draft model from a template.
"""
__author__ = 'Carolin Brune'
################################################################################
# requirements
################################################################################

import cobra
import numpy as np
import os.path
import pandas as pd
import time
import warnings

from pathlib import Path
from typing import Literal,Union

from refinegems.classes.medium import load_medium_from_db, medium_to_model
from refinegems.utility.io import load_model
from refinegems.utility.entities import resolve_compartment_names
from refinegems.curation.biomass import test_biomass_presence
from refinegems.analysis.investigate import run_memote

from refinegems.analysis.growth import MIN_GROWTH_THRESHOLD

# further required programs:
#        - MEMOTE, tested with version 0.13.0

################################################################################
# functions
################################################################################

def pid_filter(data: pd.DataFrame, pid: float) -> pd.DataFrame:
    """Filter the data based on PID threshold. Entries above the given value are kept.

    Args:
        - data (pd.DataFrame): The data from bidirectional_blast.py containing at least a 'PID' column.
        - pid (float): PID threshold value, given in percentage e.g. 80.0.

    Returns:
        pd.DataFrame: The filtered data.
    """

    data['homolog_found'] = np.where(data['PID'] > pid, 1, 0)
    counts = data['homolog_found'].value_counts()
    print(F'\t{counts[1]} set as 1 (= homologs), {counts[0]} set as 0 ')

    return data

# @TODO
def edit_template_identifiers(data:pd.DataFrame, edit:Literal['no','dot-to-underscore']) -> pd.DataFrame:
    """Edit the subject IDs to fit the gene IDs of the template model.
    Requires further extention, if needed edits are not included.

    Args:
        - data (pd.DataFrame): The data frame containing the bidirectional blastp best hits information.
        - edit (Literal['no','dot-to-underscore']): Type of edit to perform.
            Currently possible options: no, dot-to-underscore.

    Returns:
        pd.DataFrame: The (un)edited DataFrame.
    """
    

    match edit:
        case 'no':
            return data
        case 'dot-to-underscore':
            data['subject_ID'] = [x.replace('.','_') for x in data['subject_ID']]
            return data
        case _:
            mes = 'Unknown option for parameter edit. Nothing will be edited.'
            warnings.warn(mes)
            return data
        # ...........
        # @TODO
        # extendable
        # ...........


def remove_absent_genes(model:cobra.Model, genes:list[str]) -> cobra.Model:
    """Remove a list of genes from a given model.
    Note: genes that are not found in the model are skipped.

    Args:
        - model (cobra.Model): A template model to delete genes from. 
            A copy will be created before deleting.
        - genes (list[str]): Gene identifiers of genes that should be deleted.

    Returns:
        cobra.Model: A new model with the given genes deleted, if found in the original model.
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
                if variant.objective.value > MIN_GROWTH_THRESHOLD:
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


def rename_found_homologs(draft:cobra.Model, bbh:pd.DataFrame) -> cobra.Model:
    """Rename the genes in the model correnspondingly to the homologous ones found in the query.

    Args:
        - draft (cobra.Model): The draft model with the to-be-renamed genes.
        - bbh (pd.DataFrame): The table from the bidirectional_blast.py script containing the bidirectional blastp best hits information

    Returns:
        cobra.Model: The draft model with renamed genes.
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


def check_unchanged(draft:cobra.Model, bbh:pd.DataFrame) -> cobra.Model:
    """Check the genes names (more correctly, the IDs) for still existing original col_names.
    Depending on the case, decide if to keep or remove them.

    Args:
        - draft (cobra.Model): The draft model currently in the making.
        - bbh (pd.DataFrame): The table from the bidirectional_blast.py script containing the bidirectional blastp best hits information.

    Returns:
        cobra.Model: The model after the check and possible removal of genes.
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
            res = model.optimize()
            if res.objective_value > MIN_GROWTH_THRESHOLD:
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


def gen_draft_model(model:cobra.Model, bbh:pd.DataFrame, 
                    name:str, dir:str, edit:Literal['no','dot-to-underscore'], 
                    medium:str='default', namespace:Literal['BiGG']='BiGG') -> cobra.Model:
    """Generate a draft model from a template model and the results of a bidirectional blastp (blast best hits) table
    and save it as a new model.

    Args:
        - model (cobra.Model): The template model.
        - bbh (pd.DataFrame): The bidirectional blastp best hits table.
        - name (str): Name of the newly generated model.
        - dir (str): Path to the directory to save the new model in.
        - edit (Literal['no','dot-to-underscore'):  Type of edit to perform.
            Currently possible options: no, dot-to-underscore.
        - medium (str, optional):  Name of the to be loaded from the refineGEMs database or 'default' = the one
            from the template model. If given the keyword 'exchanges', will use all exchange reactions in the model as a medium.
            Defaults to 'default'.
        - namespace (Literal['BiGG'], optional): Namespace of the model. 
            Defaults to 'BiGG'.

    Returns:
        cobra.Model: The generated draft model.
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
    draft = resolve_compartment_names(draft)

    # for each object, save a note that it was added during draft construction
    for r in draft.reactions:
        r.notes['creation'] = 'via template'
    for g in draft.genes:
        g.notes['creation'] = 'via template'
    for m in draft.metabolites:
        m.notes['creation'] = 'via template'

    # save draft model
    draft.id = name
    cobra.io.write_sbml_model(draft, Path(dir,name+'_draft.xml'))

    return draft


def run(template:str, bpbbh:str, dir:str, 
        edit_names:Literal['no','dot-to-underscore']='no', 
        pid:float=80.0, name:Union[str,None]=None, 
        medium:str='default', namespace:str='BiGG', memote:bool=False):
    """Generate a draft model from a blastp best hits tsv file and a template model.

    Args:
        - template (str): Path to the file containing the template model.
        - bpbbh (str): Path to the blastp bidirectional best hits.
        - dir (str): Path to output directory.
        - edit_names (Literal['no','dot-to-underscore', optional):  Type of edit to perform.
            Currently possible options: no, dot-to-underscore. 
            Defaults to 'no'.
        - pid (float, optional): Threshold value for determining, if a gene is counted as present or absent. 
            Given in percentage, e.g. 80.0 = 80%.
            Defaults to 80.0.
        - name (Union[str,None], optional): Name of the output model. 
            If not given, takes name from filename. 
            Defaults to None.
        - medium (str, optional):  Name of the to be loaded from the refineGEMs database or 'default' = the one
            from the template model. If given the keyword 'exchanges', will use all exchange reactions in the model as a medium.
            Defaults to 'default'.
        - namespace (str, optional): Namespace of the model. 
            Defaults to 'BiGG'.
        - memote (bool, optional): Option to run memote after creating the draft model. 
            Defaults to False.

    Raises:
        ValueError: 'Edit_names value {edit_names} not in list of allowed values: no, dot-to-underscore'
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
    template_model = load_model(template,'cobra')
    # if possible, set growth function as model objective
    growth_objfunc = test_biomass_presence(template_model)
    if len(growth_objfunc) == 1:
        template_model.objective = growth_objfunc[0]
    elif len(growth_objfunc) > 1:
        mes = f'Multiple BOF detected. Chosing the following: {growth_objfunc[0]}'
        warnings.warn(mes, category=UserWarning)
        template_model.objective = growth_objfunc[0]
    else:
        mes = f'No BOF detected. Can lead to problems downstream the SPECIMEN pipeline.'
        warnings.warn(mes, category=UserWarning)

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
    draft = gen_draft_model(template_model, bbh_data, name, dir, edit_names, medium=medium, namespace=namespace)
    end = time.time()
    print(F'\ttotal time: {end - start}s')

    # -------------------
    # analyse with MEMOTE
    # -------------------

    if memote:
        memote_path = Path(dir,'step1-extension',name+'.html')
        run_memote(draft, 'html', return_res=False, save_res=memote_path, verbose=True)

    total_time_e = time.time()
    print(F'total runtime: {total_time_e-total_time_s}')

