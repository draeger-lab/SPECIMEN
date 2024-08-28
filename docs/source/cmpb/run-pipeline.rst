Run the ``CMPB`` Workflow
=========================

This page explains how to run the complete ``CMPB`` (CarveMe + ModelPolisher based) worfklow 
and how to collect the neccessary data.

For more information about the steps of the worfklow, 
see :ref:`cmpb-overview`.

``CMPB``: Quickstart 
--------------------

.. warning::

    Currently, the workflow can only be run with an already generated model as input.
    The CarveMe connection will be added in a future update.

The worfklow can either be run directly from the command line or its functions can be called from inside a Python script.
The input in both cases is a :doc:`configuration file <cmpb-config>` that contains all information needed (data file paths and parameters) to run it.

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
Missing entries will be reported when starting the worfklow.

To run the worfklow using the configuration file, use

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

The worfklow has two obligatory parameters:

- Path to a model 

    - If no model is given, the `protein_fasta` needs to be provided. The format needs to be the same as the files provides by NCBI under `<GenBank assembly>` -> `ftp` -> `<name>_translated_cds.faa.gz`

- A media configuration (from refineGEMs) for testing the model's growth

Further data can be added as available and/or needed (all are completely optional):

- The generated draft model e.g. using ``CarveMe``
- The reference sequence GFF file (for gap analysis via KEGG required, optional for CarveMe polishing)
    - Some of the gap-filling options (BioCyc, Gene) also require a GFF file, but since the type of GFF influcences the results, the input is separated from the first GFF.
- If available, the KEGG organism ID (for gap analysis via KEGG required, optional for CarveMe polishing)
- The protein FASTA of your input genome (required for lab\_strain=True, otherwise optional)
- Additional files for filling gaps: 

    - For KEGG see bullet points above 
    - Gap-filling with BioCyc requires two BioCyc SmartTables, one for the genes and one for the reactions of the organism.
    - | The gap-filling via genes uses a SwissProt database file and mapping (for more information about the setup, see ``refinegems.utility.setup.download_url``).
      | Additionally, if checking protein accession numbers against NCBI should be enabled, an email address needs to be given.

- To enable adjusting the biomass objective function using ``BOFdat``, the following information is required
    
    - Path to a file containing the full genome sequenece of your organism
    - The DNA weight fraction of your organism (experimentally determined or retrieved using literature research)
    - The enzyme/ion weight fraction of your organism (experimentally determined or retrieved using literature research)
