"""Further annotate a model.
"""
__author__ = 'Carolin Brune'
################################################################################
# requirements
################################################################################

from Bio.KEGG import REST, Enzyme
import cobra
from libsbml import *
from pathlib import Path
import subprocess
import time
from tqdm import tqdm
import urllib.error
import pandas as pd

# refinegems
from refinegems.io import load_model_cobra

# from SBOannotator import *
from SBOannotator import sbo_annotator
from ... import util

# further required programs:
#        - SBOannotator
#        - MEMOTE,  tested with version 0.13.0

################################################################################
# functions
################################################################################

def kegg_reaction_to_kegg_pathway(model, viaEC=False, viaRC=False):
    """Retrieve the KEGG pathways for existing KEGG reaction IDs, if
    they have yet to be added. Depending on the given options, only the
    reactions is searched or additional searches are started using the
    EC number and reactions class if the first search was unsuccesful.
    
    @NOTE: EC number and reaction class cover a broader set of pathways and might
    not be entirely accurate for the given reaction.

    :param model: The model to be annotated with KEGG pathway.
    :type model: cobra.Model
    :param viaEC: Option to search for KEGG pathway ID using the EC number if
        previous searches were unsuccesful. Default is False.
    :type viaEC: bool
    :param viaRC: Option to search for KEGG pathway ID using the reaction class if
        previous searches were unsuccesful. Default is False.
    :type viaEC: bool
    """

    # identify reaction with KEGG reaction annotation
    # but no KEGG pathway
    for reac in tqdm(model.reactions):
        if 'kegg.reaction' in reac.annotation and 'kegg.pathway' not in reac.annotation:

            pathways = []
            reaction = None

            # via reaction
            try:
                reaction = util.io.kegg_reaction_parser(reac.annotation['kegg.reaction'])
                if 'db' in reaction and 'kegg.pathway' in reaction['db']:
                    pathways = reaction['db']['kegg.pathway']
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


def run(model, dir, kegg_viaEC=False, kegg_viaRC=False, memote=False):
    """Further annotate a given model.

    Currently add annotations for:
    - SBO using SBOannotator

    :param model: Path to the model (sbml)
    :type model: string
    :param dir: Path to the output directory.
    :type dir: string
    :param kegg_viaEC: Option to search for KEGG pathway ID using the EC number if
        previous searches were unsuccesful. Default is False.
    :type kegg_viaEC: bool
    :param kegg_viaRC: Option to search for KEGG pathway ID using the reaction class if
        previous searches were unsuccesful. Default is False.
    :type kegg_viaEC: bool
    :param memote: Option to run memote after the annotation on the model.
        Default is False.
    :type memote: bool, optional
    """

    print('\nrefinement step 3: annotation\n################################################################################\n')

    # -----------------------
    # create output directory
    # -----------------------

    try:
        Path(F"{dir}step3-annotation/").mkdir(parents=True, exist_ok=False)
        print(F'Creating new directory {F"{dir}step3-annotation/"}')
    except FileExistsError:
        print('Given directory already has required structure.')

    # ------------------
    # add SBO annotation
    # ------------------

    print('\n# ------------------\n# add SBO annotation\n# ------------------')

    start = time.time()

    libsbml_doc = readSBML(model)
    libsbml_model = libsbml_doc.getModel()

    # note:
    #    current implementation needs the script SBOannotator.py
    #    and the create_dbs.sql files from the SBOannotator github repo
    #    runs currently only when started from 03-refinement/
    # @TODO:
    #    SBOannotator should be downloadable as a python tool
    #    if not included download from github

    sbo_annotator(libsbml_doc, libsbml_model, 'constraint-based', True, 'create_dbs', F"{dir}step3-annotation/"+libsbml_model.getId()+'_SBOannotated.xml')

    end = time.time()
    print(F'\ttime: {end - start}s')

    # reload model
    model = load_model_cobra(F"{dir}step3-annotation/"+libsbml_model.getId()+'_SBOannotated.xml')

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

    kegg_reaction_to_kegg_pathway(model, viaEC=kegg_viaEC, viaRC=kegg_viaRC)

    end = time.time()
    print(F'\ttime: {end - start}s')

    # ----------
    # save model
    # ----------
    cobra.io.write_sbml_model(model, F'{dir}step3-annotation/{model.id}_annotated.xml')

    # ---------------------------------
    # assess model quality using memote
    # ---------------------------------

    if memote:
        print('\n# -------------------\n# analyse with MEMOTE\n# -------------------')
        start = time.time()
        name = libsbml_model.getId()
        draft_path = F'{dir}step3-annotation/{name}_SBOannotated.xml'.replace(" ", "\ ")
        memote_path = F'{dir}step3-annotation/{name}_SBOannotated.html'.replace(" ", "\ ")
        subprocess.run([F'memote report snapshot --skip test_consistency --filename {memote_path} {draft_path}'], shell=True)
        end = time.time()
        print(F'\ttotal time: {end - start}s')
