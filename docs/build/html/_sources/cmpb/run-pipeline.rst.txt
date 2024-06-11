Run the CMPB Pipeline
=====================

This page explains how to run the complete CMPB (CarveMe + ModelPolisher based) pipeline 
and how to collected the neccessary data.

For more information about the steps of the pipeline, 
see :ref:`Overview of the CMPB Pipeline`.

CMPB: Quickstart 
----------------

.. warning::

    Currently, the pipeline can only be run with an already generated model as input.
    The CarveMe connection will be added in a future update.

The pipeline can either be run directly from the command line or its functions can be called from inside a Python script.
The input in both cases is a configuration file that contains all information needed (data file paths and parameters) to run it.

The configuration can be downloaded using the command line:

.. code-block:: bash

    specimen setup config -t cmpb

To download the configuration file using Python, use:

.. code-block:: python

    import specimen
    specimen.setup.download_config(filename='./my_basic_config.yaml', type='cmpb')

After downloading the configuration file, open it with an editor and change the parameters as needed.
Missing entries will be reported when starting the pipeline.

To run the pipeline using the configuration file, use

.. code-block:: bash

    specimen cmpb run "config.yaml"

on the command line or

.. code-block:: python

    specimen.cmpb.workflow.run(config_file='config.yaml')

from inside a Python script or Jupyter Notebook with "config.yaml" being the path to your configuration file.

CMPB: Collecting Data
---------------------

The pipeline has only two absolutely required parameters:

- the path to the annotated genome file (if a model is given, should be the file used to create it)
- a media configuration (from refineGEMs) for testing the model's growth

Further data can be added as available and/or needed (all are totally optional):

- the generated draft model e.g. using CarveMe
- The reference sequence GFF file (for gap analysis required, optional for CarveMe polishing)
- if available, the KEGG organism ID (for gap analysis required, optional for CarveMe polishing)
- the protein FASTA of your input genome (required for lab\_strain=True, otherwise optional)
- additional files for gapfilling

    - for KEGG see bullet points above 
    - for BioCyc, three txt files with downloaded smart tables and a protein fasta with:

         - 'Accession-2' and 'Reaction of gene' columns
         - all reaction relevant information (\*)
         - all metabolite relevant information (+)
         - protein FASTA used as input for CarveMe

    - optionally, a manually curated EXCEL sheet with information to be (potentially) added to the model

- to enable adjusting the biomass objective function using BOFdat, the following information is required
    
    - path to a file containing the full genome sequenece of your organism
    - the DNA weight fraction of your organism (experimentally determined or retrieved using literature research)
    - the enzmy/ion weight fraction of your organism (experimentally determined or retrieved using literature research)


