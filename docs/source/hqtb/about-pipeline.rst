About HQTB
==========

TODO: do not like it

The HQTB (high-quality template based) pipeline curates the model starting from the annotated 
genome of your strain of interest and an already curated, ideally very high-quality template model 
of a closely realted strain and some additional database information.

This type of pipeline aims to build upon already existing knowledge to speed up model curation
and minimize the need to perform steps again that have already been done in a similar concepts.

Overview of the HQTB Pipeline
-----------------------------

The following image shows an overrview of the steps of the pipeline:

.. image:: ../images/hqtb_pipeline-overview.png

The pipeline consists of five main steps:

.. toctree::
    :maxdepth: 3
    :numbered:

    Bidirectional BLAST <step-desc/bidirect_blast.rst>
    Draft Model Generation <step-desc/gen_draft.rst>
    Model Refinement <step-desc/refinement.rst>
    Validation <step-desc/validation.rst>
    Analysis <step-desc/analysis.rst>

The wrapper function allows the curation of multiple model sequentially with the same 
boudary conditions.

.. hint::
    Many of the steps of the pipeline can be fine tuned and turned off/on. 
    Check the configuration file for a full list of all parameters.

.. note::

    All steps of the pipeline can be run separatly via the command line or 
    the Python integration (see :ref:`Run the HTQB Pipeline`).

    All accessable function are listed in the :ref:`Contents of SPECIMEN` section.

.. footbibliography::