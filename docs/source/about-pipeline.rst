About SPECIMEN
==============

``SPECIMEN`` is a Python package for strain-specific metabolic modelling. 
Currently, its main features consist of:

- a pipeline for automated model curation based on a template and additional data
- run parts of the pipeline separatly.
- run the pipeline on a folder of different input genomes using the same parameters.

The ``SPECIMEN`` is based on the
`refineGEMS <https://github.com/draeger-lab/refinegems>`__ :footcite:p:`bauerle2023genome` toolbox.

Overview of the Pipeline
------------------------

The core of ``SPECIMEN`` is an automated pipeline, which curates a new model from an input genome, 
a high-quality template model and optional additional information.

The pipeline consists of five main steps:

.. toctree::
    :maxdepth: 3
    :numbered:

    Bidirectional BLAST <step-desc/bidirect_blast.rst>
    Draft Model Generation <step-desc/gen_draft.rst>
    Model Refinement <step-desc/refinement.rst>
    Validation <step-desc/validation.rst>
    Analysis <step-desc/analysis.rst>

.. hint::
    Many of the steps of the pipeline can be fine tuned and turned off/on. 
    Check the configuration file for a full list of all parameters.

.. note::

    All steps of the pipeline can be run separatly via the command line or 
    the Python integration (see :ref:`Run the pipeline`).

    All accessable function are listed in the :ref:`Contents of SPECIMEN` section.

.. footbibliography::