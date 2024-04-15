"""Steps of the workflow to create a GEM based on a template model.
"""

__author__ = 'Carolin Brune'

################################################################################
# requirements
################################################################################

import contextlib
import os.path
import tempfile
import warnings
import yaml

from pathlib import Path

from . import core, util

# further required programs:
#        - DIAMOND, tested with version 0.9.14 (sensitivity options fail), >2.0.4 (everything works)
#        - MEMOTE,  tested with version 0.13.0

################################################################################
# functions
################################################################################

# @TODO
# @TEST
#  -> improvements regarding new SPECIMEN update and OS independency
def run_complete(config_file:str = 'test_config.yaml'):
    """Run the complete workflow for creating a strain-specific model.

    Args:
        - config_file (str, optional): Path to the config file. 
            Defaults to 'test_config.yaml'.

    Raises:
        ValueError: Unkown file extension {extension} for file {config["subject"]["annotated_genome"]}.
    """

    # read in the configuration file
    config = util.set_up.validate_config(config_file)

    # step 0: generate output folder(s) (+ for log files)
    # current variables:
    #     config['out']['dir']: directory for output, e.g. 10-run/
    try:
        Path(config['out']['dir'],"logs").mkdir(parents=True, exist_ok=False)
        print('Creating new directory ' + Path(config["out"]["dir"],"logs"))
    except FileExistsError:
        warnings.warn('Given directory already exists. High possibility of files being overwritten.')


    # step 1: bidirectional blast
    # ---------------------------
    print('step 1/5: bidirectional blast\n\tprogress in logs->log_01_bidirectional_blast.txt')
    with open(Path(config["out"]["dir"], "logs", "log_01_bidirectional_blast.txt"),'w') as log:
        with contextlib.redirect_stdout(log):
            core.bidirectional_blast.run(config['template']['annotated_genome'],
                                         config['subject']['annotated_genome'],
                                         Path(config["out"]["dir"], "01_bidirectional_blast"),
                                         template_name=config['parameters']['bidirectional_blast']['template_name'],
                                         input_name=config['parameters']['bidirectional_blast']['input_name'],
                                         temp_header=config['parameters']['bidirectional_blast']['temp_header'],
                                         in_header=config['parameters']['bidirectional_blast']['in_header'],
                                         threads=config['performance']['threads'],
                                         sensitivity=config['parameters']['bidirectional_blast']['sensitivity'])

    # step 2: generate draft model
    # ----------------------------
    print('step 2/5: generate draft model\n\tprogress in logs->log_02_generate_draft_model.txt')
    with open(Path(config["out"]["dir"], 'logs', 'log_02_generate_draft_model.txt'),'w') as log:
        with contextlib.redirect_stdout(log):
            bpbbh = Path(config["out"]["dir"],'01_bidirectional_blast', os.path.splitext(os.path.basename(config["subject"]["annotated_genome"]))[0] + '_' + os.path.splitext(os.path.basename(config["template"]["annotated_genome"]))[0] + '_bbh.tsv')
            core.generate_draft_model.run(config['template']['model'],
                                          bpbbh,
                                          Path(config["out"]["dir"],'02_generate_draft_model'),
                                          edit_names=config['parameters']['generate_draft_model']['edit_names'],
                                          pid=config['parameters']['generate_draft_model']['pid'],
                                          name=config['out']['name'],
                                          medium=config['parameters']['generate_draft_model']['medium'],
                                          namespace=config['template']['namespace'],
                                          memote=config['out']['memote'])


    # step 3: refinement
    # ------------------
    print('step 3/5: refinement\n\tprogress in logs->log_03_refinement.txt')
    with open(Path(config["out"]["dir"],'logs','log_03_refinement.txt'),'w') as log:

        with contextlib.redirect_stdout(log):

            extension = os.path.splitext(os.path.basename(config['subject']['annotated_genome']))[1]
            match extension:
                case '.gbff':
                    fasta = Path(config["out"]["dir"],'01_bidirectional_blast','FASTA',os.path.splitext(os.path.basename(config["subject"]["annotated_genome"]))[0] + '_prot.fa')
                case '.faa':
                    fasta = config['subject']['annotated_genome']
                case _:
                    raise ValueError(F'Unkown file extension {extension} for file {config["subject"]["annotated_genome"]}.')

            core.refinement.extension.run(Path(config["out"]["dir"],'02_generate_draft_model',config["out"]["name"]+'_draft.xml'),
                                          Path(config["out"]["dir"],'01_bidirectional_blast',os.path.splitext(os.path.basename(config["subject"]["annotated_genome"]))[0]+'_info.csv'),
                                          fasta,
                                          config['data']['diamond'],
                                          Path(config["out"]["dir"]+'03_refinement'),
                                          config['data']['mnx_chem_prop'],
                                          config['data']['mnx_chem_xref'],
                                          config['data']['mnx_reac_prop'],
                                          config['data']['mnx_reac_xref'],
                                          config['data']['ncbi_map'],
                                          config['data']['ncbi_dat'],
                                          id=config['parameters']['refinement_extension']['id'],
                                          sensitivity=config['parameters']['refinement_extension']['sensitivity'],
                                          coverage=config['parameters']['refinement_extension']['coverage'],
                                          pid=config['parameters']['refinement_extension']['pid'],
                                          threads=config['performance']['threads'],
                                          exclude_dna=config['parameters']['refinement_extension']['exclude_dna'],
                                          exclude_rna=config['parameters']['refinement_extension']['exclude_rna'],
                                          memote=config['out']['memote'])
            if config['data']['universal']:
                universal = config['data']['universal']
            elif config['data']['pan-core']:
                universal = config['data']['pan-core']
            else:
                universal = None
            core.refinement.cleanup.run(Path(config["out"]["dir"],'03_refinement','step1-extension',config["out"]["name"]+'_extended.xml'),
                                        Path(config["out"]["dir"],'03_refinement'),
                                        biocyc_db=config['data']['biocyc'],
                                        check_dupl_reac = config['parameters']['refinement_cleanup']['check_dupl_reac'],
                                        check_dupl_meta = config['parameters']['refinement_cleanup']['check_dupl_meta'],
                                        remove_unused_meta = config['parameters']['refinement_cleanup']['remove_unused_meta'],
                                        remove_dupl_reac = config['parameters']['refinement_cleanup']['remove_dupl_reac'],
                                        remove_dupl_meta = config['parameters']['refinement_cleanup']['remove_dupl_meta'],
                                        universal = universal,
                                        media_path = config['parameters']['refinement_cleanup']['media_gap'],
                                        namespace = config['template']['namespace'],
                                        iterations=config['performance']['gapfilling']['iterations'],
                                        chunk_size=config['performance']['gapfilling']['chunk_size'],
                                        growth_threshold = config['parameters']['refinement_cleanup']['growth_threshold'],
                                        memote = config['out']['memote'])
            core.refinement.annotation.run(Path(config["out"]["dir"],'03_refinement','step2-clean-up',config["out"]["name"]+'_clean.xml'),
                                           Path(config["out"]["dir"],'03_refinement'),
                                           kegg_viaEC=config['parameters']['refinement_annotation']['viaEC'],
                                           kegg_viaRC=config['parameters']['refinement_annotation']['viaRC'],
                                           memote=config['out']['memote'])
            core.refinement.smoothing.run(config['subject']['full_sequence'],
                                          Path(config["out"]["dir"],'03_refinement','step3-annotation',config["out"]["name"]+'_annotated.xml'),
                                          Path(config["out"]["dir"],'03_refinement'),
                                          mcc=config['parameters']['refinement_smoothing']['mcc'],
                                          egc_solver = config['parameters']['refinement_smoothing']['egc'],
                                          namespace = config['template']['namespace'],
                                          dna_weight_frac=config['parameters']['refinement_smoothing']['dna_weight_frac'],
                                          ion_weight_frac=config['parameters']['refinement_smoothing']['ion_weight_frac'],
                                          memote=config['out']['memote'])

    # step 4: validation
    # ------------------
    print('step 4/5: validation\n\tprogress in logs->log_04_validation.txt')
    with open(Path(config["out"]["dir"],'logs','log_04_validation.txt'),'w') as log:
        with contextlib.redirect_stdout(log):
            core.validation.run(dir=config["out"]["dir"],
                                model_path =  Path(config["out"]["dir"],'03_refinement','step4-smoothing',config["out"]["name"]+'_smooth.xml'),
                                tests=None,
                                run_all=True)

    # step 5: analysis
    # ----------------
    print('step 5/5: analysis\n\tprogress in logs->log_05_analysis.txt')
    with open(Path(config["out"]["dir"],'logs','log_05_analysis.txt'),'w') as log:
        with contextlib.redirect_stdout(log):
            core.analysis.run(model_path = Path(config["out"]["dir"],'03_refinement','step4-smoothing',config["out"]["name"]+'_smooth.xml'),
                              dir = config["out"]["dir"],
                              pc_model_path = config['data']['pan-core'],
                              pc_based_on = config['parameters']['analysis']['pc_based_on'],
                              namespace=config['template']['namespace'],
                              media_path=config['parameters']['analysis']['media_analysis'],
                              test_aa_auxotrophies = config['parameters']['analysis']['test_aa_auxotrophies'],
                              pathway=config['parameters']['analysis']['pathway'])


def wrapper_pipeline(config_file:str, parent_dir:str=""):
    """Run the pipeline multiple times on a folder containing subfolders with
    subject annotated genomes and full genome sequences using the same configuration.

    Args:
        - config_file (str): config file containing the general information to run the pipeline.
            The information for the subject and output place/name can be empty.
        - parent_dir (str, optional): Path to a directory to search for subfolders containing the data. 
            Defaults to "".

    Raises:
        - ValueError: No or multiple annotated genome files found: subfolder
        - ValueError:  No or multiple full genome files found: subfolder
    """

    # load config
    with open(config_file, "r") as cfg:
        config = yaml.load(cfg, Loader=yaml.loader.FullLoader)
    # load / get subfolders to be parsed
    subfolders = [_.path for _ in os.scandir(parent_dir) if _.is_dir()]

    for subfolder in subfolders:

        current_config = config.copy()
        # check and retrieve path for annotated genome
        
        possible_anno = Path(subfolder).glob('*.gbff') + Path(subfolder).glob('*.faa')
        if len(possible_anno) == 1:
            current_anno = possible_anno[0]
        else:
            raise ValueError(F'No or multiple annotated genome files found: {subfolder}')

        # check and retrieve path for full genome
        possible_full = Path(subfolder).glob('*.fasta') + Path(subfolder).glob('*.fa') + Path(subfolder).glob('*.fna')
        if len(possible_full) == 1:
            current_full = possible_full[0]
        else:
            raise ValueError(F'No or multiple full genome files found: {subfolder}')

        # change config according to current run
        current_config['subject']['annotated_genome'] = current_anno
        current_config['subject']['full_sequence'] = current_full
        current_config['out']['dir'] = subfolder
        current_config['out']['name'] = Path(current_anno).stem

        # run pipeline with new config
        with tempfile.NamedTemporaryFile(suffix='.yaml') as temp_config:

            with open(temp_config.name, 'w') as config_stream:
                yaml.dump(current_config, config_stream)
            current_config = util.set_up.validate_config(temp_config.name)

            run_complete(temp_config)
