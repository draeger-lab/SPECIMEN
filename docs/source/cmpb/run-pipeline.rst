Run the ``CMPB`` Pipeline
=========================

This page explains how to run the complete ``CMPB`` (CarveMe + ModelPolisher based) pipeline 
and how to collect the neccessary data.

For more information about the steps of the pipeline, 
see :ref:`Overview of the CMPB Pipeline`.

``CMPB``: Quickstart 
--------------------

.. warning::

    Currently, the pipeline can only be run with an already generated model as input.
    The CarveMe connection will be added in a future update.

The pipeline can either be run directly from the command line or its functions can be called from inside a Python script.
The input in both cases is a `configuration file <cmpb-config.html>`__ that contains all information needed (data file paths and parameters) to run it.

The configuration can be downloaded using the command line:

.. code-block:: bash
    :class: copyable

    specimen setup config -t cmpb

To download the configuration file using Python, use:

.. code-block:: python
    :class: copyable

    import specimen
    specimen.setup.download_config(filename='./my_basic_config.yaml', type='cmpb')

After downloading the configuration file, open it with an editor and change the parameters as needed.
Missing entries will be reported when starting the pipeline.

To run the pipeline using the configuration file, use

.. code-block:: bash
    :class: copyable

    specimen cmpb run "config.yaml"

on the command line or

.. code-block:: python
    :class: copyable

    specimen.cmpb.workflow.run(config_file='config.yaml')

from inside a Python script or Jupyter Notebook with "config.yaml" being the path to your configuration file.

``CMPB``: Collecting Data
-------------------------

The pipeline has two obligatory parameters:

- The path to the annotated genome file (if a model is given, should be the file used to create it)
- A media configuration (from refineGEMs) for testing the model's growth

Further data can be added as available and/or needed (all are completely optional):

- The generated draft model e.g. using CarveMe
- The reference sequence GFF file (for gap analysis via KEGG required, optional for CarveMe polishing)
- If available, the KEGG organism ID (for gap analysis via KEGG required, optional for CarveMe polishing)
- The protein FASTA of your input genome (required for lab\_strain=True, otherwise optional)
- Additional files for filling gaps: 

    - For KEGG see bullet points above 
    - For BioCyc, three txt files from downloaded BioCyc SmartTables and a protein FASTA with:

         - 'Accession-2' and 'Reaction of gene' columns
         - All reaction relevant information [#]_
         - All metabolite relevant information [#]_
         - Protein FASTA used as input for CarveMe

    - Optionally, a manually curated EXCEL sheet with information to be (potentially) added to the model

- To enable adjusting the biomass objective function using BOFdat, the following information is required
    
    - Path to a file containing the full genome sequenece of your organism
    - The DNA weight fraction of your organism (experimentally determined or retrieved using literature research)
    - The enzyme/ion weight fraction of your organism (experimentally determined or retrieved using literature research)

.. [#] 'Reaction' 'Reactants of reaction' 'Products of reaction' 'EC-Number' 'KEGG Reaction' 'MetaNetX' 'Reaction-Direction' 'Spontaneous?'
.. [#] 'Compound' 'Object ID' 'Chemical Formula' 'InChI-Key' 'ChEBI'
