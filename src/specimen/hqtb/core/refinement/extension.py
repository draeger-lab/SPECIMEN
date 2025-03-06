"""Part one of the third step of the pipeline: refinement - extension.

Extends the model by mapping missing genes to a user-defined DIAMOND database.
"""

__author__ = 'Carolin Brune'

################################################################################
# requirements
################################################################################

from pathlib import Path

from tqdm import tqdm
from tqdm.auto import tqdm
tqdm.pandas()

from typing import Literal,Union

from refinegems.utility.io import load_model, write_model_to_file
from refinegems.utility.connections import run_memote


from refinegems.classes.gapfill import GeneGapFiller

# further required programs:
#        - DIAMOND, tested with version 0.9.14 (works only for certain sensitivity mode)
#                   tested with version 2.1.8 (works for all sensitivity modes for that version)

################################################################################
# functions
################################################################################

# @TEST
# @TODO : connections, e.g. input / output
# @DISCUSSION move into another module?
def run(draft:str, gff:str, fasta:str, 
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
    
    # @DISCUSSION : add more print-out / logging about which step is currently running?
    

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

    try:
        Path(dir,"manual_curation").mkdir(parents=True, exist_ok=False)
        print(F'Creating new directory {str(Path(dir,"manual_curation"))}')
    except FileExistsError:
        print('Given directory already has required structure.')

    # ----------------------
    # identify missing genes
    # ----------------------
    
    # start = time.time() # @DISCUSSION maybe add this functionality to param / report of the gapfiller
    
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
        
        # save model
        write_model_to_file(extended_model, Path(dir,'step1-extension',name+'.xml'))

        # ---------------------------------
        # assess model quality using memote
        # ---------------------------------

        if memote:
            memote_path = str(Path(dir,'step1-extension',name+'.html'))
            run_memote(draft, 'html', return_res=False, save_res=memote_path, verbose=True)
    else:
        print("\tNo missing genes found. The first refinement step, extension, will be skipped. The model will still be re-saved under the new name.")
        # save model
        write_model_to_file(draft_libsbml, Path(dir,'step1-extension',name+'.xml'))
