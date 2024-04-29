The specimen.cmd_access submodule
=================================

.. automodule:: specimen.cmd_access
   :members:
   :undoc-members:
   :show-inheritance:


After a successfull installation, ``SPECIMEN`` can be accessed via the command line
from inside the Python environment it was installed in using:

.. code-block:: bash

    specimen [OPTIONS] COMMAND [ARGS]

The following commands are available:

- ``run`` : Run the whole pipeline or a specific step.
- ``setup`` : Setup structure, data and more.


General Options
---------------

- ``--help``: Call the help page of the command.

specimen setup 
--------------

.. code:: bash

   specimen setup config 

Download a configuration file, either for the pipeline or for media.

Options:

- ``--filename/-f``: Name/Path to save the config under.
- ``--type/-t``: Type of config to download. Can be media or basic/advanced for the pipeline config.

.. code:: bash

   specimen data structure

Setup a directory with the basic structure for the data needed for the pipeline.

Options:

- ``--dir/-d``: Name/Path of the directory
- ``--chunk-size/-s``: Parameter for doenloading files from the web.

specimen run 
------------

.. code:: bash

   specimen run pipeline [CONFIG]

Run the complete pipeline with a configuration file as input.

.. code:: bash 

   specimen run wrapper [CONFIG]
   
Run the pipeline using a config on a directory containing multiple input genomes.

Options:

- ``--dir/-d``: Name/Path of the directory that contains the input.

.. code:: bash 

   specimen run bdb [TEMPLATE] [INPUT]

Run step 1: bidirectional BLAST of the pipeline. Requires the input and template genome as input.

Options:

- ``--template-name`` : Name of the annotated genome of the template, if it should not be extracted from the filename.
- ``--input-name`` : Name of the annotated genome of the input, if it should not be extracted from the filename.
- ``--temp-header``: Feature qualifier of the gbff of the template to use as header for the FASTA files. Defaults to ``protein_id``.
- ``--in-header``: Feature qualifier of the gbff of the input to use as header for the FASTA files. Defaults to ``locus_tag``.
- ``--dir/-d``: Name/Path of the directory to save the output to.
- ``--threads/-t``: Number of threads to use for DIAMOND.
- ``--sensitivity/-s``: Sensitvity mode to use for DIAMOND.

.. code:: bash

   specimen run draft [TEMPLATE] [BPBBH]

Run step 2: generate draft model of the pipeline. Requires the results of the bidirectional BLAST 
and the template model as input.

Options:

- ``--dir/-d``: Name/Path of the directory to save the output to.
- ``--edit-names``: Choose an option to change the IDs inside the bpbbh-file to fit the model, if neccessary.
- ``--pid``: Threshold value for the PID. Default to 80.0 (80 percent)
- ``--name``: Name of the output model.
- ``--medium``: Medium for the new model. Can be a name of a medium in the refineGEMS database, they keyword *default*, which uses the medium from the template or the keyword *exchanges*, which constructs a medium from all available exchange reactions.
- ``namespace``: Namespace to use for the model.
- ``--memote``: Run Memote after contructing the draft model.

.. code:: bash 

   specimen run validation [MODEL]

Run step 4: validation on a model.

Options:

- ``--dir/-d``: Name/Path of the directory to save the output to.
- ``--run-test/-t``: Specify validation tests to be run (multiple can be set). If the keyword *all* is given, runs all available tests.

.. code:: bash 

   specimen run analysis [MODEL]

Run step 5: analysis on a model.

Options:

- ``--dir``: Path to a directory for the output.
- ``--pan-core-comparison/--pcc``: Option on which feature the comparison of pan-core model and model should should be based on. Default is "id".
- ``--pan-core-model/--pcm``: Path to a pan-core model.
- ``--namespace/-n``: Namespace used by the given model. Defaults to BiGG.
- ``--media-path/--mp``: Path to a media config file. Enables growth analysis if given.
- ``--test-aa-auxotrophies/--taa``: Option to test media/model for auxotrophies.
- ``--pathway/--pathway-analysis``: Option to perform a pathway analysis using KEGG pathway identifiers.

specimen run refinement
^^^^^^^^^^^^^^^^^^^^^^^

Run the different parts of the step 3: refinement of the pipeline.

.. code:: bash

   specimen run refinement extension 

Run the first part, extension.

Required options:

- ``--draft``: Path to the draft model.
- ``--gene-list/-g``: Path to a csv file containing information on all the genes found in the annotated genome.
- ``--fasta/-f``: Path to the (protein) FASTA file containing the CDS sequences
- ``--db/--database``: Path to the database used for running DIAMOND.
- ``--mnx-chem-prop``: Path to the MetaNetX chem_prop namespace file.
- ``--mnx-chem-xref``: Path to the MetaNetX chem_xref namespace file.
- ``--mnx-reac-prop``: Path to the MetaNetX reac_prop namespace file.
- ``--mnx-reac-xref``: Path to the MetaNetX reac_xref namespace file.

Further options:

- ``--ncbi_map``: Path to the ncbi information mapping file. Optional, but recommended.
- ``--ncbi_dat``: Path to the ncbi database information file. Optional, but recommended.
- ``--dir/-d``: Path to the directory for the output (directories)
- ``--id/-i``: Name of the column of the csv file that contains the entries that were used as gene identifiers in the draft model.
- ``--sensitivity/-s``: Sensitivity mode for DIAMOND blastp run. Default is sensitive.
- ``--coverage/-c``: Threshold value for the query coverage for DIAMOND. Default is 80.0.
- ``--pid``: PID (percentage identity value) to filter the blast hist by. Default is 90.0, only hits equal or above the given value are kept.
- ``--threads/-t``: Number of threads to be used.
- ``--include_dna``: Include reactions with DNA in their name when added (developer information: True == excluded).
- ``--include_rna``: Include reactions with RNA in their name when added (developer information: True == excluded).
- ``--memote``: Use memote on the extended model.

.. code:: bash

   specimen run refinement cleanup [MODEL]

Based on a draft model, run the second part of refinement, cleanup.

Options:

- ``--dir/-d``: Path to the directory for the output (directories)
- ``--biocyc_db``: Path to the BioCyc (MetaCyc) database information file (for reactions). Optional, but recommended. Necessary for checking directionality
- ``--check_dupl_reac/--cdr``: 'Check for duplicate reactions.
- ``--check_dupl_meta/--cdm``: default='default``: Check for duplicate metabolites. Can "default" (starting point MetaNetX), exhaustive (iterate over all annotations as starting points) or "skip".
- ``--objective_function``: '--of``: Name, ID of the objective function of the model. Default is "Growth".
- ``--remove_dupl_meta/--rdm``: 'Option for removing/replacing duplicate metabolites.
- ``--remove_unused_meta/--rum``: 'Option for removing unused metabolites from the model. Only used when cdm is not skipped.
- ``--remove_dupl_reac/--rdr``: 'Option for removing duplicate reaction from the model.
- ``--universal/-u``: Path to a universal model containing reactions used for gapfilling.
- ``--media-path/--mp``: Path to a media config to use for gapfilling.
- ``--namespace/--nsp``: Namespace to use for the model.
- ``--growth_threshold/-gt``: Threshold value for a model to be considered growing.
- ``--iterations/-i``: Number of iterations for the gapfilling. If 0 is passed, uses full set of reactions instead of heuristic.
- ``--chunk_size``: Number of reactions to be tested simultaniously if using the heuristic version of gapfilling. If this is 0, heuristic will not be applied.

.. code:: bash

   specimen run refinement annotation [MODEL]
   
Run the thrid part of the refinement, annotation, on a given model.

Options:

- ``--dir``: Path to a directory for the output.
- ``--kegg-via-ec/--via-ec``: Try to map EC numbers to KEGG pathway, if KEGG reaction cannot be mapped directly.
- ``--kegg-via-rc/--via-rc``: Try to map RC numbers to KEGG pathway, if KEGG reaction cannot be mapped directly.
- ``--memote``: Use memote on the extended model.

.. code:: bash 

   specimen run refinement smoothing [MODEL]

Required Options:

- ``--genome/-g``: Path to the genome FASTA (e.g. .fna) file of your genome.

Further Options:

- ``--dir/-d``: Path to a directory for the output.
- ``--egc-solver/--egc``: String sets the type of solver to use to solve EGCs. Otherwise just reports existing EGCs.
- ``--namespace/--nsp``: Namespace of the model.
- ``--mcc``: Option to perform MassChargeCuration on the model. Can be used directly on model or as extra information. Choices are "apply","extra" and "skip". Deafult is "skip".
- ``--dna_weight_frac``: DNA macromolecular weight fraction for your organism. Default is 0.023 for Klebsiella based on Liao et al.
- ``--ion_weight_frac``: weight fraction for the coenzymes and ions. Default is 0.05 based on the default of BOFdat.
- ``--memote``: Use memote on the extended model.
