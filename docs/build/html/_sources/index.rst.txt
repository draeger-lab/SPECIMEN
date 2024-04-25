.. SPECIMEN documentation master file, created by
   sphinx-quickstart on Sun Aug  6 09:58:51 2023.

Welcome to SPECIMEN!
====================================

``SPECIMEN`` is a Python package that provides functionalities for mostly automated 
strain-specific curation of metabolic models using a high-quality template model.

Overview
------------

At core of ``SPECIMEN`` is an automated pipeline that curates a new GEM based on a genome, a template model and some additional information.
The pipeline can be accessed from the command line or from inside a Python script. An overview of the different steps of the pipeline can be seen below.

.. image:: ./images/pipeline-overview.png

Additionally, ``SPECIMEN`` allows the use of the functions / steps of the pipeline separatly.

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   installation
   about-pipeline
   run-pipeline
   specimen
   help
   dev-notes

