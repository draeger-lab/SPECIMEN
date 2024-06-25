.. SPECIMEN documentation master file, created by
   sphinx-quickstart on Sun Aug  6 09:58:51 2023.

Welcome to SPECIMEN!
====================

``SPECIMEN`` is a Python package that contains a growing collection of workflows for the 
automated curation of high-quality, ideally strain-specific, genome-scale metabolic models (GEMs).

These workflows are mainly based on the `refineGEMs <https://github.com/draeger-lab/refinegems>`__ :footcite:p:`bauerle2023genome` toolbox.

Additionally, ``SPECIMEN`` allows the use of most of the functions and 
steps of the different pipelines separatly.

Available Pipelines
-------------------

Currently, ``SPECIMEN`` includes the following pipelines (for a summary refer to :ref:`Overview of the pipelines`):

.. image:: images/buttons/cmpb.png
  :height: 0px
  :width: 0px

.. raw:: html

  <a class="reference external image-reference" href="cmpb/about-pipeline.html">
    <img src='_images/cmpb.png' alt='CMPB' title='CMPB' style="width: 30%;">
  </a>

.. image:: images/buttons/hqtb.png
  :height: 0px
  :width: 0px

.. raw:: html

  <a class="reference external image-reference" href="hqtb/about-pipeline.html">
    <img src='_images/hqtb.png' alt='HQTB' title='HQTB' style="width: 30%;">
  </a>

.. image:: images/buttons/PGAB_constr.png
  :height: 0px
  :width: 0px

.. raw:: html

  <a class="reference external image-reference" href="pipeline_idea.html">
    <img src='_images/PGAB_constr.png' alt='PGAB' title='PGAB' style="width: 30%;">
  </a>


How to Cite
-----------

.. warning::

  Coming soon


.. toctree::
   :maxdepth: 2
   :caption: Content

   installation
   overview-pipes 
   specimen
   help
   dev-notes

.. footbibliography::
