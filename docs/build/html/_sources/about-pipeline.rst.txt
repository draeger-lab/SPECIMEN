About the SPECIMEN
==================

**SPECIMEN** is a Python package for strain-specific metabolic modelling. 
Currently, its main features consist of:

- a pipeline for automated model curation based on a template and additional data
- run parts of the pipeline separatly.
- run the pipeline on a folder of different input genomes using the same parameters.

The functionalities of **SPECIMEN** are based on **refineGEMS** (link).

About the Pipeline
------------------

The core of **SPECIMEN** is an automated pipeline, which curates a new model from an input genome, a high-quality template model and optional additional information.
The pipeline consists of five steps:

1. bidirectional BLAST: perform a bidirectional BLAST using DIAMOND on the input and template genome
2. draft model generation: generate a draft model from the template model by removing or renaming genes and reactions based on the bidirectional BLAST results
3. model refinement: further refine the model

    a. extension: map not-yet-added genes to reactions and add them to the model
    b. clean-up: improve quality by checking reaction direction, duplicates, gapfilling etc.
    c. annotation: extend model annotation using e.g. SBOannotator and KEGG pathway annotations
    d. smoothing: adjust FBA parameters using MCC, BOFdat etc.

4. validation: validate the model using COBRApy validation (COBRA and libsbml validation)
5. analysis: analyse the curated model. Current options include:

    a. statistical analysis
    b. pathway analysis
    c. growth analysis, including *in silico* auxotrophy tests
    d. pan-core analysis, if a pan-core model is available

Many of the steps of the pipeline can be fine tuned and turned off/on. Check the configuration file for a full list of all parameters.

.. note::

    All steps of the pipeline can be run separatly via the command line or the Python integration.
    All accessable function are listed in the Contents of **SPECIMEN** section.
