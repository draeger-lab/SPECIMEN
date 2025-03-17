"""Further annotate a model.
"""
__author__ = 'Carolin Brune'

################################################################################
# requirements
################################################################################

from Bio.KEGG import REST, Enzyme
from libsbml import *
from pathlib import Path
from tqdm import tqdm

import cobra
import pandas as pd
import time
import urllib.error

# refinegems
from refinegems.utility.io import load_model, write_model_to_file
from refinegems.utility.connections import run_memote, run_SBOannotator
from refinegems.utility.db_access import kegg_reaction_parser
from refinegems.curation.polish import polish_annotations

################################################################################
# functions
################################################################################
# @TODO merge with refineGEMs.curation.pathways.kegg_pathways
def kegg_reaction_to_kegg_pathway(model:cobra.Model, viaEC:bool=False, viaRC:bool=False):
    """Retrieve the KEGG pathways for existing KEGG reaction IDs, if
    they have yet to be added. Depending on the given options, only the
    reactions is searched or additional searches are started using the
    EC number and reactions class if the first search was unsuccesful.

    Args:
        - model (cobra.Model): 
            The model - loaded with COBRApy - to be annotated.
        - viaEC (bool, optional): 
            Option to search for KEGG pathway ID 
            using the EC number if previous searches were unsuccesful. 
            Defaults to False.
        - viaRC (bool, optional): 
            Option to search for KEGG pathway ID 
            using the reaction class if previous searches were unsuccesful. 
            Defaults to False.
    """

    # identify reaction with KEGG reaction annotation
    # but no KEGG pathway
    for reac in tqdm(model.reactions):
        if 'kegg.reaction' in reac.annotation and 'kegg.pathway' not in reac.annotation:

            pathways = []
            reaction = None

            # via reaction
            try:
                if isinstance(reac.annotation['kegg.reaction'],list):
                    for annotation in reac.annotation['kegg.reaction']:
                        reaction = kegg_reaction_parser(annotation)
                        if reaction is not None and 'db' in reaction and 'kegg.pathway' in reaction['db']:
                            if isinstance(reaction['db']['kegg.pathway'],list):
                                pathways.extend(reaction['db']['kegg.pathway'])
                            else:
                                pathways.append(reaction['db']['kegg.pathway'])
                else:
                    reaction = kegg_reaction_parser(reac.annotation['kegg.reaction'])
                    if reaction is not None and 'db' in reaction and 'kegg.pathway' in reaction['db']:
                        if isinstance(reaction['db']['kegg.pathway'],list):
                                pathways.extend(reaction['db']['kegg.pathway'])
                        else:
                            pathways.append(reaction['db']['kegg.pathway'])
            except urllib.error.HTTPError:
                print(F'HTTPError: {reac.id}, {reac.annotation["kegg.reaction"]}')
            except ConnectionResetError:
                print(F'ConnectionResetError: {reac.id}, {reac.annotation["kegg.reaction"]}')
            except urllib.error.URLError:
                print(F'URLError: {reac.id}, {reac.annotation["kegg.reaction"]}')

            # via reaction class
            # can lead to some additional classes, as the RC are not as strictly defined as
            # the reactions themselves
            if viaRC and len(pathways) == 0 and not pd.isnull(reaction) and 'rc' in reaction:
                try:
                    collect = False
                    for kegg_rc in reaction['rc']:
                        kegg_rc = REST.kegg_get(kegg_rc)
                        kegg_rc = kegg_rc.read()
                        for line in kegg_rc.split('\n'):
                            if line:
                                if line.startswith('PATHWAY'):
                                    collect = True
                                    pathways.append(line.replace('PATHWAY','',1).strip().split(' ')[0])
                                elif collect == True and line[0] != '/':
                                    if line[0].isupper():
                                        collect = False
                                    else:
                                        pathways.append(line.strip().split(' ')[0])
                except urllib.error.HTTPError:
                    print(F'HTTPError: {reac.id}, {reaction["rc"]}')
                except ConnectionResetError:
                    print(F'ConnectionResetError: {reac.id}, {reaction["rc"]}')
                except urllib.error.URLError:
                    print(F'URLError: {reac.id}, {reaction["rc"]}')

            # via EC
            # seems really sketchy to do it this way, as ONE EC number
            # can include MANY reactions for different pathways
            if viaEC and len(pathways) == 0:
                if 'ec-code' in reac.annotation:
                    ec_code = reac.annotation['ec-code']
                elif pd.isnull(reaction) and 'db' in reaction and 'ec-code' in reaction['db']:
                    ec_code = reaction['db']['ec-code']
                else:
                    ec_code = '-'
                if ec_code != '-' and isinstance(ec_code,str):
                    try:
                        kegg_ec = REST.kegg_get(F'ec:{ec_code}')
                        kegg_ec = Enzyme.read(kegg_ec)
                        if len(kegg_ec.pathway) == 0 or kegg_ec.pathway == None:
                            pass
                        else:
                            for i in kegg_ec.pathway:
                                pathways.append(i[1])
                    except urllib.error.HTTPError:
                        print(F'HTTPError: {reac.id}, {ec_code}')
                    except ConnectionResetError:
                        print(F'ConnectionResetError: {reac.id}, {ec_code}')
                    except urllib.error.URLError:
                        print(F'URLError: {reac.id}, {ec_code}')
                else:
                    print(F'No EC number: {reac.id}')

            # add pathway annotation to reaction if found
            if len(pathways) != 0:
                reac.annotation['kegg.pathway'] = pathways

# @TEST
def run(model:str, dir:str, kegg_viaEC:bool=False, 
        kegg_viaRC:bool=False, memote:bool=False):
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
    

    print('\nrefinement step 3: annotation\n################################################################################\n')

    # -----------------------
    # create output directory
    # -----------------------

    try:
        Path(dir,"step3-annotation").mkdir(parents=True, exist_ok=False)
        print(F'Creating new directory {str(Path(dir,"step3-annotation"))}')
    except FileExistsError:
        print('Given directory already has required structure.')

    # --------------
    # Load the model
    # --------------

    model = load_model(model, 'libsbml')

    # ------------------
    # Polish annotations
    # ------------------

    print('\n# ----------------------------------\n# polish annotations\n# ----------------------------------')
    start = time.time()

    model = polish_annotations(model, True, str(Path(dir,'step3-annotation',model.getId()+'_annotations_polished.xml')))

    end = time.time()
    print(F'\ttime: {end - start}s')

    # ------------------
    # add SBO annotation
    # ------------------

    print('\n# ------------------\n# add SBO annotation\n# ------------------')

    start = time.time()

    model = run_SBOannotator(model)
    write_model_to_file(model, str(Path(dir,'step3-annotation',model.getId()+'_SBOannotated.xml')))

    end = time.time()
    print(F'\ttime: {end - start}s')

    # reload model
    model = load_model(str(Path(dir,'step3-annotation',model.getId()+'_SBOannotated.xml')), 'cobra')

    # ................................................................
    # @EXTENDABLE
    #    possible to add more functions using model as in- and output
    #    to further annotated the model
    # ................................................................

    # ----------------
    # add KEGG pathway
    # ----------------
    print('\n# ----------------\n# add KEGG pathway\n# ----------------')

    start = time.time()

    # @TODO Compare to pathways.kegg_pathways in refineGEMs
    kegg_reaction_to_kegg_pathway(model, viaEC=kegg_viaEC, viaRC=kegg_viaRC)

    end = time.time()
    print(F'\ttime: {end - start}s')

    # ----------
    # save model
    # ----------
    cobra.io.write_sbml_model(model, Path(dir,'step3-annotation',model.id+'_annotated.xml'))

    # ---------------------------------
    # assess model quality using memote
    # ---------------------------------

    if memote:
        start = time.time()
        name = model.id
        memote_path = str(Path(dir,'step3-annotation',name+'_annotated.html'))
        run_memote(model, 'html', return_res=False, save_res=memote_path, verbose=True)
        end = time.time()
        print(F'\ttotal time: {end - start}s')