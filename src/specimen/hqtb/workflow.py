"""Functions to run the workflow to create a GEM based on a high-quality template model.
"""

__author__ = "Carolin Brune"

################################################################################
# requirements
################################################################################

import contextlib
import os.path
import tempfile
import warnings
import yaml

from datetime import date
from pathlib import Path

from refinegems.utility.io import mimic_genbank

from . import core
from .. import util

# further required programs:
#        - DIAMOND, tested with version 0.9.14 (sensitivity options fail), >2.0.4 (everything works)
#        - MEMOTE,  tested with version >=0.13.0

################################################################################
# functions
################################################################################


def run(config_file: str = "test_config.yaml"):
    """Run the complete workflow for creating a strain-specific model.

    Args:
        - config_file (str, optional):
            Path to the config file.
            Defaults to 'test_config.yaml'.

    Raises:
        - ValueError: Unkown file extension.
    """

    # read in the configuration file
    config = util.set_up.validate_config(config_file)

    # step 0: generate output folder(s) (+ for log files)
    # current variables:
    #     config['general']['dir']: directory for output, e.g. 10-run/
    try:
        Path(config["general"]["dir"], "logs").mkdir(parents=True, exist_ok=False)
        print("Creating new directory " + str(Path(config["general"]["dir"], "logs")))
    except FileExistsError:
        warnings.warn(
            "Given directory already exists. High possibility of files being overwritten."
        )

    # set up model name
    # -----------------
    if (
        config["general"]["authorinitials"] is not None
        and config["general"]["organism"] is not None
        and config["general"]["strainid"] is not None
    ):
        modelname = (
            "i"
            + config["general"]["organism"]
            + str(config["general"]["strainid"])
            + config["general"]["authorinitials"]
            + str(date.today().year).removeprefix("20")
        )
    elif config["general"]["modelname"] is not None:
        modelname = config["general"]["modelname"]
    else:
        print(
            "No values given for the standard name for a model. Default name will be used."
        )
        modelname = "model_" + str(date.today().year).removeprefix("20")

    # step 1: bidirectional blast
    # ---------------------------
    print(
        "step 1/5: bidirectional blast\n\tprogress in logs->log_01_bidirectional_blast.txt"
    )
    with open(
        Path(config["general"]["dir"], "logs", "log_01_bidirectional_blast.txt"), "w"
    ) as log:
        with contextlib.redirect_stdout(log):
            core.bidirectional_blast.run(
                config["template"]["annotated_genome"],
                config["subject"]["annotated_genome"],
                Path(config["general"]["dir"], "01_bidirectional_blast"),
                template_name=config["parameters"]["bidirectional_blast"][
                    "template_name"
                ],
                input_name=config["parameters"]["bidirectional_blast"]["input_name"],
                temp_header=config["parameters"]["bidirectional_blast"]["temp_header"],
                in_header=config["parameters"]["bidirectional_blast"]["in_header"],
                threads=config["performance"]["threads"],
                sensitivity=config["parameters"]["bidirectional_blast"]["sensitivity"],
            )

    # step 2: generate draft model
    # ----------------------------
    print(
        "step 2/5: generate draft model\n\tprogress in logs->log_02_generate_draft_model.txt"
    )
    with open(
        Path(config["general"]["dir"], "logs", "log_02_generate_draft_model.txt"), "w"
    ) as log:
        with contextlib.redirect_stdout(log):
            bpbbh = Path(
                config["general"]["dir"],
                "01_bidirectional_blast",
                os.path.splitext(
                    os.path.basename(config["subject"]["annotated_genome"])
                )[0]
                + "_"
                + os.path.splitext(
                    os.path.basename(config["template"]["annotated_genome"])
                )[0]
                + "_bbh.tsv",
            )
            core.generate_draft_model.run(
                config["template"]["model"],
                bpbbh,
                Path(config["general"]["dir"], "02_generate_draft_model"),
                edit_names=config["parameters"]["generate_draft_model"]["edit_names"],
                pid=config["parameters"]["generate_draft_model"]["pid"],
                name=modelname,
                medium=config["parameters"]["generate_draft_model"]["medium"],
                namespace=config["template"]["namespace"],
                memote=config["general"]["memote"],
            )

    # step 3: refinement
    # ------------------
    print("step 3/5: refinement\n\tprogress in logs->log_03_refinement.txt")
    with open(
        Path(config["general"]["dir"], "logs", "log_03_refinement.txt"), "w"
    ) as log:

        with contextlib.redirect_stdout(log):

            # step 3.1: extension
            # ...................
            
            # create a GenBank format FASTA for the extension step (needed for the GapFiller)
            fasta_path = mimic_genbank(config["subject"]["annotated_genome"], config["subject"]["gff"],
                                       str(Path(config["general"]["dir"]))) # @ASK any ides for a better place for this file?
            
            core.refinement.extend(
                draft=Path(
                    config["general"]["dir"],
                    "02_generate_draft_model",
                    modelname + "_draft.xml",
                ),
                gff=config["subject"]["gff"],
                fasta=fasta_path, 
                db=config["data"]["diamond"],
                dir=Path(config["general"]["dir"] + "03_refinement"),
                ncbi_mapping=config["data"]["ncbi_map"],
                email=config["parameters"]["general"]["email"],
                sensitivity=config["parameters"]["refinement_extension"]["sensitivity"],
                coverage=config["parameters"]["refinement_extension"]["coverage"],
                pid=config["parameters"]["refinement_extension"]["pid"],
                threads=config["performance"]["threads"],
                threshold_add_reacs=config["parameters"]["refinement_extension"][
                    "threshold_add_reacs"
                ],
                prefix=config["parameters"]["general"]["idprefix"],
                namespace=config["parameters"]["general"]["namespace"],
                formula_check=config["parameters"]["refinement_extension"][
                    "formula-check"
                ],
                exclude_dna=config["parameters"]["refinement_extension"]["exclude-dna"],
                exclude_rna=config["parameters"]["refinement_extension"]["exclude-rna"],
                memote=config["general"]["memote"],
            )
            # step 3.2: cleanup
            # .................

            if config["data"]["universal"]:
                universal = config["data"]["universal"]
            elif config["data"]["pan-core"]:
                universal = config["data"]["pan-core"]
            else:
                universal = None

            if config["parameters"]["refinement_cleanup"]["GeneGapFiller"]:
                gene_gapfiller_params = {
                    "prefix": config["parameters"]["refinement_cleanup"]["idprefix"],
                    "type_db": config["parameters"]["refinement_cleanup"][
                        "GeneGapFiller parameters"
                    ]["type"],
                    "fasta": config["parameters"]["refinement_cleanup"][
                        "GeneGapFiller parameters"
                    ]["fasta"],
                    "dmnd_db": config["parameters"]["refinement_cleanup"][
                        "GeneGapFiller parameters"
                    ]["dmnd-database"],
                    "map_db": config["parameters"]["refinement_cleanup"][
                        "GeneGapFiller parameters"
                    ]["database-mapping"],
                    "mail": config["parameters"]["general"]["email"],
                    "check_NCBI": config["parameters"]["refinement_cleanup"][
                        "GeneGapFiller parameters"
                    ]["check-NCBI"],
                    "threshold_add_reacs": config["parameters"]["refinement_cleanup"][
                        "threshold_add_reacs"
                    ],
                    "sens": config["parameters"]["refinement_cleanup"][
                        "GeneGapFiller parameters"
                    ]["sensitivity"],
                    "cov": config["parameters"]["refinement_cleanup"][
                        "GeneGapFiller parameters"
                    ]["coverage"],
                    "t": config["performance"]["threads"],
                    "pid": config["parameters"]["refinement_cleanup"][
                        "GeneGapFiller parameters"
                    ]["percentage identity"],
                    "formula_check": config["parameters"]["refinement_cleanup"][
                        "formula-check"
                    ],
                    "exclude_dna": config["parameters"]["refinement_cleanup"][
                        "exclude-dna"
                    ],
                    "exclude_rna": config["parameters"]["refinement_cleanup"][
                        "exclude-rna"
                    ],
                    "gff": config["parameters"]["refinement_cleanup"][
                        "GeneGapFiller parameters"
                    ]["gff"],
                }
            else:
                gene_gapfiller_params = None

            core.refinement.cleanup(
                Path(
                    config["general"]["dir"],
                    "03_refinement",
                    "step1-extension",
                    modelname + "_extended.xml",
                ),
                Path(config["general"]["dir"], "03_refinement"),
                run_gene_gapfiller=gene_gapfiller_params,
                biocyc_db=config["data"]["biocyc"],
                check_dupl_reac=config["parameters"]["refinement_cleanup"][
                    "check_dupl_reac"
                ],
                check_dupl_meta=config["parameters"]["refinement_cleanup"][
                    "check_dupl_meta"
                ],
                remove_unused_meta=config["parameters"]["refinement_cleanup"][
                    "remove_unused_meta"
                ],
                remove_dupl_reac=config["parameters"]["refinement_cleanup"][
                    "remove_dupl_reac"
                ],
                remove_dupl_meta=config["parameters"]["refinement_cleanup"][
                    "remove_dupl_meta"
                ],
                universal=universal,
                media_path=config["parameters"]["refinement_cleanup"]["media_gap"],
                namespace=config["template"]["namespace"],
                iterations=config["performance"]["gapfilling"]["iterations"],
                chunk_size=config["performance"]["gapfilling"]["chunk_size"],
                growth_threshold=config["parameters"]["refinement_cleanup"][
                    "growth_threshold"
                ],
                memote=config["general"]["memote"],
            )

            # step 3.3: annotation
            # ....................

            core.refinement.annotate(
                Path(
                    config["general"]["dir"],
                    "03_refinement",
                    "step2-clean-up",
                    modelname + "_clean.xml",
                ),
                Path(config["general"]["dir"], "03_refinement"),
                kegg_viaEC=config["parameters"]["refinement_annotation"]["viaEC"],
                kegg_viaRC=config["parameters"]["refinement_annotation"]["viaRC"],
                memote=config["general"]["memote"],
            )

            # step 3.4: smoothing
            # ...................

            core.refinement.smooth(
                config["subject"]["full_sequence"],
                Path(
                    config["general"]["dir"],
                    "03_refinement",
                    "step3-annotation",
                    modelname + "_keggpathways.xml",
                ),
                Path(config["general"]["dir"], "03_refinement"),
                mcc=config["parameters"]["refinement_smoothing"]["mcc"],
                egc_solver=config["parameters"]["refinement_smoothing"]["egc"],
                namespace=config["template"]["namespace"],
                dna_weight_frac=config["parameters"]["refinement_smoothing"][
                    "dna_weight_frac"
                ],
                ion_weight_frac=config["parameters"]["refinement_smoothing"][
                    "ion_weight_frac"
                ],
                memote=config["general"]["memote"],
            )

    # step 4: validation
    # ------------------
    print("step 4/5: validation\n\tprogress in logs->log_04_validation.txt")
    with open(
        Path(config["general"]["dir"], "logs", "log_04_validation.txt"), "w"
    ) as log:
        with contextlib.redirect_stdout(log):
            core.validation.run(
                dir=config["general"]["dir"],
                model_path=Path(
                    config["general"]["dir"],
                    "03_refinement",
                    "step4-smoothing",
                    modelname + "_smooth.xml",
                ),
                tests=config["parameters"]["validation"]["tests"], 
                run_all=config["parameters"]["validation"]["run_all"], 
            )

    # step 5: analysis
    # ----------------
    print("step 5/5: analysis\n\tprogress in logs->log_05_analysis.txt")
    with open(
        Path(config["general"]["dir"], "logs", "log_05_analysis.txt"), "w"
    ) as log:
        with contextlib.redirect_stdout(log):
            core.analysis.run(
                model_path=Path(
                    config["general"]["dir"],
                    "03_refinement",
                    "step4-smoothing",
                    modelname + "_smooth.xml",
                ),
                dir=config["general"]["dir"],
                pc_model_path=config["data"]["pan-core"],
                pc_based_on=config["parameters"]["analysis"]["pc_based_on"],
                namespace=config["template"]["namespace"],
                media_path=config["parameters"]["analysis"]["media_analysis"],
                test_aa_auxotrophies=config["parameters"]["analysis"][
                    "test_aa_auxotrophies"
                ],
                pathway=config["parameters"]["analysis"]["pathway"],
            )


def wrapper(config_file: str, parent_dir: str = ""):
    """Run the pipeline multiple times on a folder containing subfolders with
    subject annotated genomes and full genome sequences using the same configuration.

    Args:
        - config_file (str):
            config file containing the general information to run the pipeline.
            The information for the subject and output place/name can be empty.
        - parent_dir (str, optional):
            Path to a directory to search for subfolders containing the data.
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

        possible_anno = Path(subfolder).glob("*.gbff") + Path(subfolder).glob("*.faa")
        if len(possible_anno) == 1:
            current_anno = possible_anno[0]
        else:
            raise ValueError(
                f"No or multiple annotated genome files found: {subfolder}"
            )

        # check and retrieve path for full genome
        possible_full = (
            Path(subfolder).glob("*.fasta")
            + Path(subfolder).glob("*.fa")
            + Path(subfolder).glob("*.fna")
        )
        if len(possible_full) == 1:
            current_full = possible_full[0]
        else:
            raise ValueError(f"No or multiple full genome files found: {subfolder}")

        # change config according to current run
        current_config["subject"]["annotated_genome"] = current_anno
        current_config["subject"]["full_sequence"] = current_full
        current_config["general"]["dir"] = subfolder
        current_config["general"]["modelname"] = Path(current_anno).stem

        # run pipeline with new config
        with tempfile.NamedTemporaryFile(suffix=".yaml") as temp_config:

            with open(temp_config.name, "w") as config_stream:
                yaml.dump(current_config, config_stream)
            current_config = util.set_up.validate_config(temp_config.name)

            run(temp_config)
