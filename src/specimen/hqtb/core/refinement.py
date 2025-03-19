"""Function for the refinement step (step 3) of the HQTB workflow.

Refinement includes four main parts:
- Part 1: extension
- Part 2: cleanup
- Part 3: annotation
- Part 4: smoothing 
"""

__author__ = 'Carolin Brune'

################################################################################
# requirements
################################################################################

import cobra
import os
import pandas as pd
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
from refinegems.curation.curate import resolve_duplicates
from refinegems.curation.biomass import check_normalise_biomass
from refinegems.curation.pathways import set_kegg_pathways
from refinegems.curation.polish import polish_annotations
from refinegems.utility.connections import adjust_BOF, perform_mcc, run_memote, run_SBOannotator
from refinegems.utility.io import load_model, write_model_to_file
from refinegems.utility.util import test_biomass_presence

# further required programs:
#        - DIAMOND, tested with version 0.9.14 (works only for certain sensitivity mode)
#                   tested with version 2.1.8 (works for all sensitivity modes for that version)


################################################################################
# variables
################################################################################

GGF_REQS = {'prefix', 'type_db', 'fasta', 'dmnd_db', 'map_db', 'mail', 'check_NCBI',
            'threshold_add_reacs', 'sens', 'cov', 't', 'pid', 'formula_check',
            'exclude_dna','exclude_rna','gff'}

################################################################################
# functions
################################################################################

# Part 1: extension
# -----------------

def extend(draft:str, gff:str, fasta:str, 
        db:str, dir:str, 
        # mapping to NCBI
        ncbi_mapping:Union[Path,str,None]=None,
        email:Union[None,str]=None, 
        # params for DIAMOND
        sensitivity:Literal['sensitive', 'more-sensitive', 'very-sensitive','ultra-sensitive']='more-sensitive', 
        coverage:float=95.0, 
        pid:float=90.0, 
        threads:int=2, 
        # param for adding entities to model
        threshold_add_reacs:int = 5,
        prefix:str='specimen',
        namespace:str='BiGG',
        formula_check:Literal['none','existence','wildcard','strict']='existence',
        exclude_dna:bool=True, exclude_rna:bool=True, 
        # validation
        memote:bool=False):
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
    
    if not sensitivity in ['sensitive','more-sensitive','very-sensitive','ultra-sensitive']:
        raise ValueError(F'Unknown sensitivity mode {sensitivity}. Choose from: sensitive, more-sensitive, very-sensitive, ultra-sensitive')

    print('\nrefinement step 1: extension\n################################################################################\n')

    # -----------------------
    # create output directory
    # -----------------------

    try:
        Path(dir,"step1-extension").mkdir(parents=True, exist_ok=False)
        print(F'Creating new directory {str(Path(dir,"step1-extension"))}')
    except FileExistsError:
        print('Given directory already has required structure.')

    # ----------------------
    # identify missing genes
    # ----------------------

    # load model
    draft_libsbml = load_model(draft, 'libsbml')
    draft_cobra = load_model(draft, 'cobra')
    
    name = f'{draft_cobra.id}_extended'
    
    # set up GapFiller
    gp = GeneGapFiller()

    # identify missing genes
    gp.find_missing_genes(gff, draft_libsbml) 

    if hasattr(gp,'missing_genes') and len(gp.missing_genes) >= 1:
        
        # --------------------------
        # identify missing reactions
        # --------------------------
        
        kwargs = {'sens': sensitivity, 'cov': coverage, 't': threads, 
                  'pid': pid, 'outdir': dir}
        gp.find_missing_reactions(draft_cobra, prefix, 'user', fasta, db, 
                                  ncbi_mapping, email, False, threshold_add_reacs, 
                                  **kwargs)

        # ------------
        # extend model
        # ------------
        
        kwargs = {'namespace': namespace, 'idprefix': prefix, 'formula_check': formula_check, 
                  'exclude_dna': exclude_dna, 'exclude_rna': exclude_rna} 
        extended_model = gp.fill_model(draft_libsbml, **kwargs)
        
        # save GapFiller report
        gp.report(Path(dir,"step1-extension"))
        
        # save model
        write_model_to_file(extended_model, Path(dir,'step1-extension',name+'.xml'))

        # ---------------------------------
        # assess model quality using memote
        # ---------------------------------

        if memote:
            print('\nRunning memote ...\n------------------\n')
            memote_path = str(Path(dir,'step1-extension',name+'.html'))
            run_memote(draft, 'html', return_res=False, save_res=memote_path, verbose=True)
    else:
        print("\tNo missing genes found. The first refinement step, extension, will be skipped. The model will still be re-saved under the new name.")
        # save model
        write_model_to_file(draft_libsbml, Path(dir,'step1-extension',name+'.xml'))


# Part 2: cleanup
# ---------------

# similar to the one in refineGEMs, but add an extra check to 
# skip reactions created using the template
def check_direction(model:cobra.Model,data_file:str) -> cobra.Model:
    """Check the direction of newly created reactions (01-extention) by searching for matching MetaCyc,
    KEGG and MetaNetX IDs as well as EC number in a downloaded BioCyc (MetaCyc)
    database table (need to contain at least the following columns:
    Reactions (MetaCyc ID),EC-Number,KEGG reaction,METANETX,Reaction-Direction.

    ..note::

        Checks only creations that do not contain the notes['creation'] == 'via template',
        assuming the template was well curated.

    Args:
        - model (cobra.Model): 
            The GEM containing the reactions to be check for direction.
        - data_file (str): 
            Path to the MetabCyc (BioCyc) smart table.

    Returns:
        cobra.Model: 
            The updated model.
    """

    # create MetaCyc table
    # --------------------
    data = pd.read_csv(data_file, sep='\t')
    # rewrite the columns into a better comparable/searchable format
    data['KEGG reaction'] = data['KEGG reaction'].str.extract(r'.*>(R\d*)<.*')
    data['METANETX']      = data['METANETX'].str.extract(r'.*>(MNXR\d*)<.*')
    data['EC-Number']     = data['EC-Number'].str.extract(r'EC-(.*)')

    # check direction
    # --------------------
    for r in model.reactions:
            # entry from template, assumed to be already curated
            if 'creation' in r.notes and 'via template' == r.notes['creation']:
                continue
            # newly created entry, check direction with BioCyc
            else:
                direction = None
                # easy case: metacyc is already (corretly) annotated
                if 'metacyc.reaction' in r.annotation and len(data[data['Reactions'] == r.annotation['metacyc.reaction']]) != 0:
                    direction = data[data['Reactions'] == r.annotation['metacyc.reaction']]['Reaction-Direction'].iloc[0]
                    r.notes['BioCyc direction check'] = F'found {direction}'
                # complicated case: no metacyc annotation
                else:
                    annotations = []

                    # collect matches
                    if 'kegg.reaction' in r.annotation and r.annotation['kegg.reaction'] in data['KEGG reaction'].tolist():
                        annotations.append(data[data['KEGG reaction'] == r.annotation['kegg.reaction']]['Reactions'].tolist())
                    if 'metanetx.reaction' in r.annotation and r.annotation['metanetx.reaction'] in data['METANETX'].tolist():
                        annotations.append(data[data['METANETX'] == r.annotation['metanetx.reaction']]['Reactions'].tolist())
                    if 'ec-code' in r.annotation and r.annotation['ec-code'] in data['EC-Number'].tolist():
                        annotations.append(data[data['EC-Number'] == r.annotation['ec-code']]['Reactions'].tolist())

                    # check results
                    # no matches
                    if len(annotations) == 0:
                        r.notes['BioCyc direction check'] = 'not found'

                    # matches found
                    else:
                        # built intersection
                        intersec = set(annotations[0]).intersection(*annotations)
                        # case 1: exactly one match remains
                        if len(intersec) == 1:
                            entry = intersec.pop()
                            direction = data[data['Reactions'] == entry]['Reaction-Direction'].iloc[0]
                            r.annotation['metacyc.reaction'] = entry
                            r.notes['BioCyc direction check'] = F'found {direction}'

                        # case 2: multiple matches found -> inconclusive
                        else:
                            r.notes['BioCyc direction check'] = F'found, but inconclusive'

                # update direction if possible and needed
                if not pd.isnull(direction):
                    if 'REVERSIBLE' in direction:
                        # set reaction as reversible by setting default values for upper and lower bounds
                        r.lower_bound = -1000.
                    elif 'RIGHT-TO-LEFT' in direction:
                        # invert the default values for the boundaries
                        r.lower_bound = -1000.
                        r.upper_bound = 0.
                    else:
                        # left to right case is the standart for adding reactions
                        # = nothing left to do
                        continue
    return model


def cleanup(model:str, dir:str, 
        biocyc_db:str=None, 
        run_gene_gapfiller:Union[None,dict]=None,
        check_dupl_reac:bool = False,
        check_dupl_meta:bool = 'default',
        remove_unused_meta:bool = False, 
        remove_dupl_reac:bool = False, 
        remove_dupl_meta:bool = False,
        universal:str = None, 
        media_path:str = None, 
        namespace:Literal['BiGG']='BiGG', 
        growth_threshold:float = 0.05,
        iterations:int=3, chunk_size:int=10000,
        memote:bool = False):
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
            The dictionary needs to contain the keys saved in :py:data:`~specimen.hqtb.core.refinement.cleanup.GGF_REQS`.
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
    if not check_dupl_meta in ['default','skip','exhaustive']:
        raise ValueError('Unknown option {check_dupl_meta} for checking duplicate metabolite. Use one of: default, skip, exhaustive')
    
    if run_gene_gapfiller:
        if not GGF_REQS.issubset(run_gene_gapfiller.keys()):
            raise KeyError('At least one parameter for the GeneGapFiller is missing. Re-check your input for run_gene_gapfiller')

    # -------------
    # start program
    # -------------
    print('\nrefinement step 2: clean-up\n################################################################################\n')

    # -----------------------
    # create output directory
    # -----------------------

    try:
        Path(dir,"step2-clean-up").mkdir(parents=True, exist_ok=False)
        print(F'Creating new directory {str(Path(dir,"step2-clean-up"))}')
    except FileExistsError:
        print('Given directory already has required structure.')

    model = cobra.io.read_sbml_model(model)

    # --------------------
    # check direction
    # --------------------

    if biocyc_db:

        print('\n# --------------------\n# check direction\n# --------------------')

        start = time.time()

        # check direction
        model = check_direction(model,biocyc_db)

        end = time.time()
        print(F'\ttime: {end - start}s')

    # ----------
    # gapfilling
    # ----------

    print('\n# ----------\n# gapfilling\n# ----------')
    start = time.time()
    
    # gapfilling
    ############
        
    # GeneGapFiller 
    # -------------
    if run_gene_gapfiller:
        # load model with libsbml
        libmodel = None
        with NamedTemporaryFile(suffix='.xml', delete=False) as tmp:
            write_model_to_file(model,tmp.name)
            libmodel = load_model(tmp.name,'libsbml')
        os.remove(tmp.name)
        
        # run the gene gap filler
        ggf = GeneGapFiller()
        ggf.find_missing_genes(run_gene_gapfiller['gff'],
                               libmodel)
        ggf.find_missing_reactions(model,
                                   prefix=run_gene_gapfiller['prefix'],
                                   type_db=run_gene_gapfiller['type_db'],
                                   fasta = run_gene_gapfiller['fasta'],
                                   dmnd_db = run_gene_gapfiller['dmnd_db'],
                                   map_db = run_gene_gapfiller['map_db'],
                                   mail = run_gene_gapfiller['mail'],
                                   check_NCBI = run_gene_gapfiller['check_NCBI'],
                                   threshold_add_reacs = run_gene_gapfiller['threshold_add_reacs'],
                                   outdir = Path(dir,"step2-clean-up"),
                                   sens = run_gene_gapfiller['sens'],
                                   cov = run_gene_gapfiller['cov'],
                                   t = run_gene_gapfiller['t'],
                                   pid = run_gene_gapfiller['pid'] 
                               )
        libmodel = ggf.fill_model(libmodel, 
                                  namespace = namespace,
                                  idprefix = run_gene_gapfiller['prefix'],
                                  formula_check = run_gene_gapfiller['formula_check'],
                                  exclude_dna = run_gene_gapfiller['exclude_dna'],
                                  exclude_rna = run_gene_gapfiller['exclude_rna']
                                )
        
        write_model_to_file(libmodel,Path(dir,"step2-clean-up",'after_2ndGF_nogff.xml'))
        
        # re-load model with cobrapy
        with NamedTemporaryFile(suffix='.xml', delete=False) as tmp:
            write_model_to_file(libmodel,tmp.name)
            model = load_model(tmp.name,'cobra')
        os.remove(tmp.name)
        
        ggf.report(Path(dir,"step2-cleanup"))

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
        universal_model = load_model(universal,'cobra')
        # run gapfilling
        model = multiple_cobra_gapfill(model,universal_model,media_list,namespace,iterations=iterations, chunk_size=chunk_size, growth_threshold=growth_threshold)

    end = time.time()
    print(F'\ttime: {end - start}s')

    # -----------------
    # resove duplicates
    # -----------------
    print('\n# -----------------\n# resolve duplicates\n# -----------------')
    start = time.time()

    model = resolve_duplicates(model,
                               check_reac = check_dupl_reac,
                               check_meta = check_dupl_meta,
                               remove_unused_meta = remove_unused_meta,
                               remove_dupl_reac = remove_dupl_reac,
                               replace_dupl_meta = remove_dupl_meta)

    end = time.time()
    print(F'\ttime: {end - start}s')

    # ---------------------
    # dead ends and orphans
    # ---------------------
    # (currently no removal of dead ends and orphans as they may be interesting
    # for manual curation)

    # ----------
    # save model
    # ----------
    name = F'{model.id}_clean'
    cobra.io.write_sbml_model(model, Path(dir,'step2-clean-up',name+'.xml'))

    # ---------------------------------
    # assess model quality using memote
    # ---------------------------------
        
    if memote:
        memote_path = str(Path(dir,'step2-clean-up',name+'.html'))
        run_memote(model, 'html', return_res=False, save_res=memote_path, verbose=True)


# Part 3: annotation
# ------------------

def annotate(model:str, dir:str, kegg_viaEC:bool=False, 
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

    set_kegg_pathways() # @TODO adjust for parameters in refinegems

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
        

# Part 4: smoothing
# -----------------

def smooth(genome:str,model:str,dir:str,mcc='skip',
        egc_solver:None|Literal['greedy']=None,
        namespace:Literal['BiGG']='BiGG',
        dna_weight_frac=0.023,ion_weight_frac=0.05, 
        memote=False):
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
    
    print('\nrefinement step 4: smoothing\n################################################################################\n')

    # -----------------------
    # create output directory
    # -----------------------

    try:
        Path(dir,"step4-smoothing").mkdir(parents=True, exist_ok=False)
        print(F'Creating new directory {str(Path(dir,"step4-smoothing"))}')
    except FileExistsError:
        print('Given directory already has required structure.')

    try:
        Path(dir,"manual_curation").mkdir(parents=True, exist_ok=False)
        print(F'Creating new directory {str(Path(dir,"manual_curation"))}')
    except FileExistsError:
        print('Given directory already has required structure.')

    # ---------
    # load data
    # ---------
    model = cobra.io.read_sbml_model(model)

    # ---------------
    # mass and charge
    # ---------------

    if mcc == 'apply':
        print('\n# ----------------------------------\n# mass and charge curation (applied)\n# ----------------------------------')
        start = time.time()
        model = perform_mcc(model,Path(dir,"step4-smoothing"))
        end = time.time()
        print(F'\ttime: {end - start}s')
    elif mcc == 'extra':
        print('\n# --------------------------------\n# mass and charge curation (extra)\n# --------------------------------')
        start = time.time()
        model = perform_mcc(model,Path(dir,"manual_curation"),False)
        end = time.time()
        print(F'\ttime: {end - start}s')
    elif mcc == 'skip':
        print('\n# ------------------------\n# mass and charge curation\n# ------------------------\n\tskipped')
    else:
        warnings.warn(F'Unknown option {mcc} for Mass and Charge Curation. Usage of MCC will be skipped.')

    # ----------------------------------
    # check for energy generating cycles
    # ----------------------------------
    print('\n# ---------------------------------------------\n# # check for energy generating cycles\n# ---------------------------------------------')
    start = time.time()

    match egc_solver:
        # greedy solver
        case 'greedy':
            print('Using GreedyEGCSolver...')
            solver = egcs.GreedyEGCSolver()
            results = solver.solve_egcs(model,namespace=namespace) # automatically uses c,e as compartments 
            if results:
                for k,v in results.items():
                    print(f'\t{k}: {v}')
        
        # no solver = EGCs will only be reported
        case _:
            solver = egcs.EGCSolver()
            print(f'\tFound EGCs:\n')
            print(f'\t{solver.find_egcs(model,with_reacs=True,namespace=namespace)}') # automatically uses c,e as compartments 

    end = time.time()
    print(F'\ttime: {end - start}s')

    # --------------------------
    # biomass objective function
    # --------------------------
    # adjust the BOF to the current genome

    print('\n# ----------\n# adjust BOF\n# ----------')
    start = time.time()

    with tempfile.NamedTemporaryFile(suffix='.xml', delete=False) as temp_model:
        # generate an up-to-date model xml-file
        cobra.io.write_sbml_model(model,temp_model.name)
        # update BOF
        pos_bofs = test_biomass_presence(model)
        if pos_bofs:
            model.reactions.get_by_id(pos_bofs[0]).reaction = adjust_BOF(genome, temp_model.name, model, dna_weight_frac, ion_weight_frac)
            # optimise BOF(s)
            model = check_normalise_biomass(model)
        else:
            # create new BOF
            bof_reac = Reaction('Biomass_BOFdat')
            bof_reac.name = 'Biomass objective function created by BOFdat'
            model.add_reactions([bof_reac])
            model.reactions.get_by_id(bof_reac).reaction = adjust_BOF(genome, temp_model.name, model, dna_weight_frac, ion_weight_frac)
       
            # optimise BOF(s)
            model = check_normalise_biomass(model)
    os.remove(temp_model.name)

    end = time.time()
    print(F'\ttime: {end - start}s')

    # ----------------
    # save final model
    # ----------------
    print('\n# ----------\n# save model\n# ----------')
    model_name = F'{model.id}_smooth'
    outname = Path(dir,'step4-smoothing',model_name+".xml")
    print(F'\tsaving to: {outname}')
    cobra.io.write_sbml_model(model,outname)

    # ---------------------------------
    # assess model quality using memote
    # ---------------------------------

    if memote:
        memote_path = str(Path(dir,'step4-smoothing',model_name+'.html'))
        run_memote(model, 'html', return_res=False, save_res=memote_path, verbose=True)

