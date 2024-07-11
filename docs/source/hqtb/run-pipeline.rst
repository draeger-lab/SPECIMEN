Run the ``HQTB`` Pipeline
=========================

This page explains how to run the complete ``HTQB`` (high-quality template based) pipeline 
and how to collected the neccessary data.

For more information about the steps of the pipeline, 
see :ref:`overview-hqtb`.

``HQTB``: Quickstart
--------------------

The pipeline can either be run directly from the command line or its functions can be called from inside a Python script.
The input in both cases is a configuration file that contains all information needed (data file paths and parameters) to run it.

The `configuration <hqtb-config.html>`__ can be downloaded using the command line:

.. code-block:: bash
    :class: copyable

    specimen setup config

This downloads a basic version with minimal parameters suitable for beginners. To download the advanced version that allows to adjust more parameters,
add :code:`-t hqtb-advanced` to the command. For further options, refer to the manual page (:code:`specimen setup config --help`).

To download the configuration file using Python, use:

.. code-block:: python
    :class: copyable

    import specimen
    specimen.util.set_up.download_config(filename='./my_basic_config.yaml', type='hqtb-basic')

As with the command line access, the type can be changed to ``hqtb-advanced``.

After downloading the configuration file, open it with an editor and change the parameters as needed.
Missing entries will be reported when starting the pipeline.

To run the pipeline using the configuration file, use

.. code-block:: bash
    :class: copyable

    specimen hqtb run pipeline "config.yaml"

on the command line or

.. code-block:: python
    :class: copyable

    specimen.hqtb.workflow.run_complete(config_file='config.yaml')

from inside a Python script or Jupyter Notebook with "config.yaml" being the path to your configuration file.

.. note::

    Additionally, the pipeline can be run with a wrapper to susequently build multiple models for different genomes using the same parameters.
    The wrapper can be accessed using :code:`specimen hqtb run wrapper "config.yaml"` or :code:`specimen.workflow.wrapper_pipeline(config_file='/User/path/to/config.yaml', parent_dir="./")`.


``HQTB``: Collecting Data
-------------------------

If you are just starting a new project and do not have all the data ready to go, you can use the setup function of
``SPECIMEN`` to help you collect the data you need.

.. code-block:: python
    :class: copyable

    specimen.util.set_up.build_data_directories('your_folder_name')

| The function above creates the following directory structure for your project.
| The 'contains' column lists what is supposed to be inside the according folder. 
  The tags manual/semi/automated report how these files are added to the folder (automated = by the setup function; semi = multiple steps neccessary, some by the program, some by the user; manual = by the user).
  The tags required/optional report whether this input is necessary to run the pipeline or if it is an optional input.

.. table::
    :align: center 

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

.. note::

    Regarding the annotated_genomes folder, the program currently only supports 
    the file types ``GBFF`` and ``FAA`` + ``FNA`` (from the NCBI and PROKKA annotation pipelines respectively)
    as genome annotation formats.

Further details for collecting the data:

- `BioCyc <https://biocyc.org/>`__:

    - Downloading a smart table from BioCyc requires a subscription.
    - The SmartTable needs to have the columns 'Reactions', 'EC-Number', 'KEGG reaction', 'METANETX' and 'Reaction-Direction'.

- RefSeq

    - One way to build a DIAMOND reference database is to download a set of reference sequences from the NCBI database, e.g. in the **FAA** format.
    - Use the function :code:`specimen.util.util.create_DIAMOND_db_from_folder('/User/path/input/directory', '/User/Path/for/output/', name = 'database', extention = 'faa')` to create a DIAMOND database
    - To speed up the mapping, create an additional mapping file from the e.g. ``GBFF`` files from NCBI using :code:`specimen.util.util.create_NCBIinfo_mapping('/User/path/input/directory', '/User/Path/for/output/', extention = 'gbff')`
    - To ensure correct mapping to KEGG, an additional information file can be created by constructing a CSV file with the following columns: 'NCBI genome', 'organism', 'locus_tag' (only the part until the seperator '_', the part, that is the same for all locus tags) and 'KEGG.organism'

        - The information of the first three columns can be taken from the previous two steps while
        - For the last column the user needs to check, if the genomes have been entered into KEGG and have an organism identifier.
        - This file is purely optional for running the pipeline but potentially leads to better results.

- medium:   

    The media, either for analysis or gap filling can be entered into the pipeline via a config file. 
    The same media file can be used for both or one file for each step can be entered into the pipeline. 
    The config files are from the `refineGEMs <https://github.com/draeger-lab/refinegems/tree/dev-2>`__ :footcite:p:`bauerle2023genome` toolbox and access its in-build medium database. 
    Additionally, the config files allow for manual adjustment / external input.

    An examplary config file can be accessed using the following command:

    .. code-block:: python
        :class: copyable

        download_config(filename='my_media_config.yaml', type='media')

    Or via the command line (additional name can be added using the flag :code:`-f <name>`):

    .. code-block:: bash
        :class: copyable
        
        specimen setup config -t media

.. note::
    The setup can be done via the command line as well, refer to :code:`specimen setup --help`.

.. footbibliography::