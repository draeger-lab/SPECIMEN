About HQTB
==========

The high-quality template based (``HTQB``) workflow curates the model starting from the annotated 
genome of your strain of interest and an already curated, ideally very high-quality template model 
of a closely related strain (species) and additional database information.

This type of workflow aims to build upon already existing knowledge to speed up model curation
and minimize the need to perform steps again that have already been done in a similar context.

.. _overview-hqtb:

Overview of the ``HQTB`` Workflow
---------------------------------

The following image shows an overrview of the steps of the workflow:

.. image:: ../images/hqtb_pipeline-overview.png

The workflow consists of five main steps:

.. toctree::
    :maxdepth: 3

    Step 1: Bidirectional BLAST <step-desc/bidirect_blast.rst>
    Step 2: Draft Model Generation <step-desc/gen_draft.rst>

- | Step 3: Model Refinement 

  | Refine the previously generated draft model to make the model more complete and strain-specific. In other words, fitting it more closely to the input genome.

    - :doc:`Part 1: Extension <step-desc/refine-parts/extension>`
    - :doc:`Part 2: Clean-up <step-desc/refine-parts/cleanup>`
    - :doc:`Part 3: Annotation <step-desc/refine-parts/annot>`
    - :doc:`Part 4: Smoothing <step-desc/refine-parts/smooth>`

.. toctree::
    :maxdepth: 3

    Step 4: Validation <step-desc/validation.rst>
    Step 5: Analysis <step-desc/analysis.rst>


.. toctree::
    :hidden:
    :maxdepth: 2

    Part 1: Extension <step-desc/refine-parts/extension.rst>
    Part 2: Clean-up <step-desc/refine-parts/cleanup.rst>
    Part 3: Annotation <step-desc/refine-parts/annot.rst>
    Part 4: Smoothing <step-desc/refine-parts/smooth.rst>



The wrapper function allows the curation of multiple models sequentially using the same 
boudary parameters.

.. hint::
    Many of the steps of the workflow can be fine tuned and turned off/on. 
    Check the :doc:`configuration file <hqtb-config>` for a full list of all parameters.

.. note::

    All steps of the workflow can be run separatly via the command line or 
    the Python integration (see :doc:`run-pipeline`).

    All accessible functions are listed in the :ref:`Contents of SPECIMEN` section.

.. footbibliography::