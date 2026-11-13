The specimen.cmd_access submodule
=================================

.. automodule:: specimen.cmd_access
   :members:
   :undoc-members:
   :show-inheritance:

After a successful installation, ``SPECIMEN`` can be accessed via the command line
from inside the Python environment it was installed in using:

.. code-block:: bash

    specimen [OPTIONS] COMMAND [ARGS]

The following commands are available:

- ``cmpb`` : Workflow for GEM curation based on CarveMe or a CarveMe model and ModelPolisher.
- ``hqtb`` : Workflow for GEM curation based on a high-quality template.
- ``setup`` : Setup structure, data and more.


General Options
---------------

- ``--help/-h``: Show top-level help. Subcommands support ``--help``.
- ``--version``: Show the installed version (top-level command only).

Uppercase arguments below are required positional arguments.

specimen setup
--------------

.. code:: bash

   specimen setup config

Download a configuration file, either for the workflow or for media.

Options:

- ``--filename/-f``: Name/Path to save the config under. Defaults to ``config.yaml``.
- ``--type/-t``: Type of config to download. Choices are ``hqtb-basic`` (default), ``hqtb-advanced``, ``hqtb-defaults``, ``media`` and ``cmpb``.

.. code:: bash

   specimen setup data-structure WORKFLOW

Setup a directory with the basic structure for the data needed for the workflow.

Argument:

- ``WORKFLOW``: The workflow: ``hqtb`` or ``cmpb``. The aliases ``"high-quality template based"`` and ``"carveme modelpolisher based"`` are also accepted.

Options:

- ``--dir/-d``: Name/Path of the directory
- ``--chunk-size/-s``: Parameter for downloading files from the web.

.. code:: bash

   specimen setup make-ncbi-mapping FOLDER

Generate a mapping file for quicker gene gap-filling.

Argument:

- ``FOLDER``: Path to the input directory to search for files for the mapping.

Options:

- ``--out/-o``: Path to the output directory.
- ``--extension/-e``: File extension of the annotated genomes to search for. Defaults to ``gbff``.

specimen hqtb
-------------

.. code:: bash

   specimen hqtb run CONFIG

Run the complete workflow with a configuration file as input.

.. code:: bash

   specimen hqtb wrapper CONFIG

Run the workflow using a config on a parent directory containing one subdirectory
per input genome, each with its annotated and full genome files.

Options:

- ``--directory/-d``: Path to the parent directory containing the input subdirectories.

.. code:: bash

   specimen hqtb bdb TEMPLATE INPUT

Run step 1: bidirectional BLAST of the workflow.

Arguments:

- ``TEMPLATE``: Path to the annotated genome of the template.
- ``INPUT``: Path to the annotated genome file of the organism-of-interest.

Options:

- ``--template-name`` : Name of the annotated genome of the template, if it should not be extracted from the filename.
- ``--input-name`` : Name of the annotated genome of the input, if it should not be extracted from the filename.
- ``--temp-header``: Feature qualifier of the gbff of the template to use as header for the FASTA files. Defaults to ``protein_id``.
- ``--in-header``: Feature qualifier of the gbff of the input to use as header for the FASTA files. Defaults to ``locus_tag``.
- ``--dir/-d``: Name/Path of the directory to save the output to.
- ``--threads/-t``: Number of threads to use for DIAMOND.
- ``--sensitivity/-s``: Sensitivity mode to use for DIAMOND.

.. code:: bash

   specimen hqtb draft TEMPLATE BPBBH

Run step 2: generate draft model of the workflow.

Arguments:

- ``TEMPLATE``: Path to the template model.
- ``BPBBH``: Path to the BLASTp bidirectional best hits results file from step 1.

Options:

- ``--dir``: Name/Path of the directory to save the output to.
- ``--edit-names``: Choose an option to change the IDs inside the bpbbh-file to fit the model, if necessary.
- ``--pid``: Threshold value for the PID. Defaults to 80.0 (80 percent)
- ``--name``: Name of the output model.
- ``--medium``: Medium for the new model. Can be a name of a medium in the refineGEMS database, the keyword *default*, which uses the medium from the template or the keyword *exchanges*, which constructs a medium from all available exchange reactions.
- ``--namespace/--nsp``: Namespace to use for the model. Only ``BiGG`` is accepted; defaults to ``BiGG``.
- ``--memote``: Run MEMOTE after constructing the draft model.

.. code:: bash

   specimen hqtb validation MODEL

Run step 4 of the ``HQTB`` workflow, validation, on a model.

Argument:

- ``MODEL``: Path to the model to be validated, e.g. output of the refinement step.

Options:

- ``--dir/-d``: Name/Path of the directory to save the output to.
- ``--run-test/-t``: Specify validation tests to be run (multiple can be set). Defaults to ``all``. Supported values are ``all``, ``cobra`` and ``modelpolisher``; ``all`` runs all available tests.

.. code:: bash

   specimen hqtb analysis MODEL

Run step 5: analysis on a model.

Arguments:

- ``MODEL``: Path the model to be analysed.

Options:

- ``--dir``: Path to a directory for the output.
- ``--pan-core-comparison/--pcc``: Option on which feature the comparison of pan-core model and model should be based on. Default is "id".
- ``--pan-core-model/--pcm``: Path to a pan-core model.
- ``--namespace/-n``: Namespace used by the given model. Only ``BiGG`` is accepted; defaults to ``BiGG``.
- ``--media-path/--mp``: Path to a media config file. Enables growth analysis if given.
- ``--test-aa-auxotrophies/--taa``: Option to test media/model for auxotrophies.
- ``--pathway/--pathway-analysis``: Option to perform a pathway analysis using KEGG pathway identifiers.

specimen hqtb refinement
^^^^^^^^^^^^^^^^^^^^^^^^

Run the different parts of the step 3: refinement of the ``HQTB`` workflow.

.. code:: bash

   specimen hqtb refinement extension

Run the first part, extension.

Required options:

- ``--draft``: Path to the draft model.
- ``--gff``: Path to the gff file of the draft model.
- ``--fasta/-f``: Path to the (protein) FASTA file containing the CDS sequences
- ``--db/--database``: Path to the database used for running DIAMOND.

Further options:

- ``--ncbi-map``: Path to the NCBI information mapping file. Optional, but recommended.
- ``--dir/-d``: Path to the directory for the output (directories)
- ``--mail/-m``: Mail address for accessing NCBI.
- ``--sensitivity/-s``: Sensitivity mode for DIAMOND BLASTp run. Default is sensitive.
- ``--coverage/-c``: Threshold value for the query coverage for DIAMOND. Default is 80.0.
- ``--pid``: PID (percentage identity value) to filter the blast hits by. Default is 95.0, only hits equal or above the given value are kept.
- ``--threads/-t``: Number of threads to be used. Defaults to 2.
- ``--threshold_add_reacs/--thres``: Maximum number of mapped reaction for an EC number allowed while adding them.
- ``--prefix``: Prefixes for entities where the ID needs to be constructed from scratch. Defaults to `HQTBext`.
- ``--namespace/-n``: Namespace to use for the model IDs. Only ``BiGG`` is accepted; this command sets no default.
- ``--formula-check``: Chemical formula acceptance level: ``none``, ``existence`` (default), ``wildcard`` or ``strict``.
- ``--include-dna``: Include reactions with DNA in their name. Excluded by default.
- ``--include-rna``: Include reactions with RNA in their name. Excluded by default.
- ``--memote``: Use MEMOTE on the extended model.

.. code:: bash

   specimen hqtb refinement cleanup MODEL DIR

Run the second part of refinement: cleanup and optional gap-filling.

Arguments:

- ``MODEL``: Path to an existing model file.
- ``DIR``: Path to an existing output directory.

Options:

- ``--reac-direc``: Path to a BioCyc SmartTable for checking reaction direction.
- ``--universal/-u``: Path to a universal model for COBRApy gap-filling.
- ``--media-path/--mp``: Path to a media configuration for gap-filling.
- ``--growth-threshold/--gt``: Growth threshold for gap-filling. Defaults to 0.05.
- ``--check-dupl-reac/--cdr``: Check for duplicate reactions.
- ``--check-dupl-meta/--cdm``: Duplicate metabolite checking mode: ``default`` (default), ``exhaustive`` or ``skip``.
- ``--remove-unused-meta/--rum``: Remove unused metabolites.
- ``--remove-dupl-reac/--rdr``: Remove duplicate reactions.
- ``--remove-dupl-meta/--rdm``: Remove duplicate metabolites.
- ``--namespace/--nsp``: Namespace of the model. Only ``BiGG`` is accepted; defaults to ``BiGG``.
- ``--iterations/-i``: Number of gap-filling iterations. Defaults to 3.
- ``--chunk-size/-c``: Chunk size for gap-filling. Defaults to 10000.
- ``--memote``: Run MEMOTE on the processed model.

For COBRApy gap-filling, provide both a universal model and a media configuration.
Boolean flags are disabled by default.

.. code:: bash

   specimen hqtb refinement annotation MODEL

Run the third part of the ``HQTB`` refinement step, annotation, on a given model.

Argument:

- ``MODEL``: Path to the model to be annotated.

Options:

- ``--dir/-d``: Path to a directory for the output.
- ``--kegg-via-ec/--via-ec/--ec``: Try to map EC numbers to KEGG pathway, if KEGG reaction cannot be mapped directly.
- ``--kegg-via-rc/--via-rc/--rc``: Try to map RC numbers to KEGG pathway, if KEGG reaction cannot be mapped directly.
- ``--memote``: Use MEMOTE on the extended model.

.. code:: bash

   specimen hqtb refinement smoothing MODEL

Argument:

- ``MODEL``: Path to the model to be smoothed.

Required Options:

- ``--genome/-g``: Path to the genome FASTA (e.g. .fna) file of your genome.

Further Options:

- ``--dir/-d``: Path to a directory for the output.
- ``--egc-solver/--egc``: Use ``greedy`` to solve energy-generating cycles (EGCs). If omitted, only report existing EGCs.
- ``--limit``: Maximum number of cores for EGC processing. Defaults to all available cores.
- ``--chunksize``: Chunk size for parallel EGC processing. Defaults to 100.
- ``--namespace/--nsp``: Namespace of the model. Only ``BiGG`` is accepted; defaults to ``BiGG``.
- ``--mcc``: Option to perform MassChargeCuration on the model. Can be used directly on model or as extra information. Choices are "apply","extra" and "skip". Default is "skip".
- ``--dna_weight_frac``: DNA macromolecular weight fraction for your organism. Default is 0.023 for Klebsiella based on Liao et al.
- ``--ion_weight_frac``: weight fraction for the coenzymes and ions. Default is 0.05 based on the default of BOFdat.
- ``--memote``: Use MEMOTE on the extended model.

specimen cmpb
-------------

.. code:: bash

   specimen cmpb run CONFIG

Run the complete CMPB workflow with a configuration file as input.
