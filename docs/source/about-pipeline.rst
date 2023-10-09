About the SPECIMEN
==================

**SPECIMEN** is a Python package for strain-specific metabolic modelling. Currently, its main features consist of:

- a pipeline for automated model curation based on a template and additional data
- a module for creating, storing and manipulating media to test a model growth on
- a set of additional functions for preparing and collecting data


About the Pipeline
------------------

The core of **SPECIMEN** is an automated pipeline, which curates a new model from an input genome, a high-quality template model and optional additional information.
The pipeline consists of five steps:

1. bidirectional BLAST: perform a bidirectional BLAST using DIAMOND on the input and template genome
2. draft model generation: generate a draft model from the template model by removing or renaming genes and reactions based on the bidirectional BLAST results
3. model refinement: further refine the model

    a. extension: map not-yet-added genes to reactions and add them to the model
    b. clean-up: improve quality by checking reaction direction, duplicates, gapfilling etc.
    c. annotation: extend model annotation using SBOannotator and KEGG pathway annotations
    d. smoothing: adjust FBA parameters using MCC, BOFdat etc.

4. validation: validate the model using COBRApy validation (COBRA and libsbml validation)
5. analysis: analyse the curated model. Current options include:

    a. statistical analysis
    b. pathway analysis
    c. growth analysis, including *in silico* auxotrophy tests
    d. pan-core analysis, if pan-core model is available

Many of the steps of the pipeline can be fine tuned and turned off/on. Check the configuration file for a full list of all parameters.

.. note::

    All steps of the pipeline can be run separatly via the command line or the Python integration.
    All accessable function are listed in the Contents of **SPECIMEN** section.


About the Media
---------------

The tool provides a small database for media to test the growth on as well as the functions to create new ones either from scratch or from the available ones.
Currently, the database includes the following media:

- M9: minimal medium
- M9_o2s: minimal medium, oxygen substituted with superoxide
- LB: lysogeny broth
- LB_o2s: lysogeny broth, oxygen substituted with superoxide
- SNM3: synthetic nasal medium
- SNM3_o2s: synthetic nasal medium, oxygen substituted with superoxide
- VMH-EUavg: VMH gut medium, EU average diet
- VMH-vegan: VMH gut medium, vegan diet
- VMH-unhealthy: VMH gut medium, unhealthy diet
- VMH-hiFloC: VMH gut medium, high fat low carb diet
