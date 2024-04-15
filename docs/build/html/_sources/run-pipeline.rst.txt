Run the Pipeline
=================

Quickstart
-----------

The pipeline can either be run directly from the command line or its functions can be called from inside a Python script.
The input in both cases is a configuration file that contains all information needed (data file paths and parameters) to run it.

The configuration can be downloaded using the command line:

.. code-block:: bash

    specimen setup config

This downloads a basic version with minimal parameters suitable for beginners. To download the advanced version that allows to adjust more parameters,
add :code:`-t advanced` to the command. For further options, refer to the manual page (:code:`specimen setup config --help`).

To download the configuration file using Python, use:

.. code-block:: python

    import specimen
    specimen.setup.download_config(filename='./my_basic_config.yaml', type='basic')

As with the command line access, the type can be changed to "advanced".

After downloading the configuration file, open it with an editor and change the parameters as needed.
Missing entries will be reported when starting the pipeline.

To run the pipeline using the configuration file, use

.. code-block:: bash

    specimen run pipeline "config.yaml"

on the command line or

.. code-block:: python

    specimen.workflow.run_complete(config_file='config.yaml')

from inside a Python script or Jupyter Notebook with "config.yaml" being the path to your configuration file.

.. note::

    Additionally, the pipeline can be run with a wrapper to susequently build multiple models for different genome using the same parameters.
    To wrapper can be accessed using :code:`specimen run wrapper "config.yaml"` or :code:`specimen.workflow.wrapper_pipeline(config_file='/User/path/to/config.yaml', parent_dir="./")`.


Collecting Data
---------------

If you are just starting a new project and do not have all the data ready to go, you can use the setup function of
**SPECIMEN** to help you collect the data you need.

.. code-block:: python

    specimen.set_up.build_data_directories('your_folder_name')

The function above creates the following directory structure for your project:

+--------------------+------------------------------+---------------------+
| folder             | contains                     | tags                |
+====================+==============================+=====================+
|| annotated_genomes || template + input            || manual, required   |
||                   || annotated + full            ||                    |
||                   || genome files                ||                    |
+--------------------+------------------------------+---------------------+
| BioCyc             | BioCyc smart table           | manual, optional    |
+--------------------+------------------------------+---------------------+
| medium             | media config, external media | manual, optional    |
+--------------------+------------------------------+---------------------+
| MetaNetX           | MetaNetX mappings            | automated, required |
+--------------------+------------------------------+---------------------+
| pan-core-models    | pan-core models              | manual, optional    |
+--------------------+------------------------------+---------------------+
|| RefSeqs           || DIAMOND database            || semi, required     |
||                   || for BLAST                   ||                    |
+--------------------+------------------------------+---------------------+
| template-models    | template models              | manual, required    |
+--------------------+------------------------------+---------------------+
| universal-models   | universal models             | manual, optional    |
+--------------------+------------------------------+---------------------+

In the contains columns it is listed what is supposed to be inside that folder.
The tags manual/semi/automated report how these are added to the folder (automated = by the setup function, manual = by the user).
The tags report/optional report whether this input is necessary to run the pipeline or if it is an optional input.

.. note::

    Regarding the annotated genomes, the program currently only supports the file types **GBFF** and **FAA** + **FNA**.

Further details for collecting the data:

- BioCyc:

    - downloading a smart table from BioCyc requires a subscription
    - the smart table needs to have the columns Reactions, EC-Number, KEGG reaction, METANETX and Reaction-Direction

- RefSeqs

    - one way to builf a DIAMOND reference database is to download a set of reference sequences from the NCBI database, e.g. in the **FAA** format
    - use the function :code:`specimen.util.util.create_DIAMOND_db_from_folder('/User/path/input/directory', '/User/Path/for/output/', name = 'database', extention = 'faa')` to create a DIAMOND database
    - to speed up the mapping, create an additional mapping file from the e.g. **GBFF** files from NCBI using :code:`specimen.util.util.create_NCBIinfo_mapping('/User/path/input/directory', '/User/Path/for/output/', extention = 'gbff')`
    - to ensure correct mapping to KEGG, an additional information file can be created by constructing a CSV file with the following columns: NCBI genome, organism, locus_tag (start) and KEGG.organism

        - the information of the first three columns can be taken from the previous two steps while
        - the last column the user needs to check, if the genomes have been entered into KEGG and have an organism identifier
        - this file is purely optional for running the pipeline but potentially leads to better results

- medium:   

    The media, either for analysis or gapfilling can be entered into the pipeline via a config file (each).
    The config files are from the **refineGEMS** toolbox and access the in-build medium database of refinegems 
    and additionally allow for manual adjustment / external input.

    A examplary config file can be accessed using the following command:

    .. code-block:: python

        download_config(filename='my_media_config.yaml', type='media')

    Or via the command line (additional name can be added using the flag :code:`-f <name>`):

    .. code-block:: bash
        
        specimen setup config -t media

.. note::
    The setup can be done via the command line as well, refer to :code:`specimen setup --help`.
