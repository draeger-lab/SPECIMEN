.. SPECIMEN documentation master file, created by
   sphinx-quickstart on Sun Aug  6 09:58:51 2023.

Welcome to SPECIMEN!
====================

``SPECIMEN`` is a Python package that contains a growing collection of workflows for the 
automated curation of high-quality, ideally strain-specific, genome-scale metabolic models (GEMs).

These workflows are mainly based on the `refineGEMS <https://github.com/draeger-lab/refinegems>`__ :footcite:p:`bauerle2023genome` toolbox.

Additionally, ``SPECIMEN`` allows the use of most of the functions and 
steps of the different pipelines separatly.

Available Pipelines
-------------------

Currently, ``SPECIMEN`` includes the following pipelines:

- | ``CMPB`` - CarveMe + ModelPolisher based: 
  | This pipeline curated the model based on a CarveMe draft model and additional input.
- | ``HQTB`` - high-quality template based:
  | Curates a model based on the annotated genome, a high-quality template model and additional database information.


.. toctree::
   :maxdepth: 2
   :caption: Content

   installation
   overview-pipes 
   specimen
   help
   dev-notes

