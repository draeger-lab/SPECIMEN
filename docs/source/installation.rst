Installation
==============

Installation via GitHub/pip
---------------------------
Download the ``SPECIMEN`` repository `here <https://github.com/cb-Hades/klebsiella-pipeline/tree/main>`_ 
and run the command :code:`pip install -e .` inside the top-level directory.

.. hint::

    To additionally install all packages necessary for working with the docs, 
    install the tool with the extra `docs`, e.g. 

    .. code:: console
        :class: copyable

        pip install -e ".[docs]" --config-settings editable_mode=strict  

For all options to run error-free, the following tools need to be installed additionally:

- `DIAMOND, version 2.0.4 or higher <https://github.com/bbuchfink/diamond>`_
- `EntrezDirect <https://www.ncbi.nlm.nih.gov/books/NBK179288/>`_, if no NCBI mapping has been created beforehand

Afterwards, ``SPECIMEN`` can either be accessed via the command line or via importing the package into a Python script.

``SPECIMEN`` can also be used via Docker. To build the docker image, firtsly clone the repository:

.. code:: console
   :class: copyable

   git clone "https://github.com/draeger-lab/specimen.git"

Then change into the directory and build the image:

.. code:: console
   :class: copyable

   cd specimen \
   docker build -t specimen .

The default command executed by the image is ``specimen -h`` and provides the help information for the CLI of 
``SPECIMEN``.

.. code:: console
   :class: copyable

   docker run specimen

To use the image interactively and open a bash shell, run the following command:

.. code:: console
   :class: copyable

   docker run -it --entrypoint bash specimen

To use the image for specific commands, you can simply use every of the CLI commands as entrypoint. 
For example, to run the CMPB pipeline, use:

.. code:: console
   :class: copyable

   docker run --name <container_name> -v <user_folder>:/sp_cont specimen cmpb run ./path/to/CMPB_config.yaml

.. hint::

    It is advised to install ``SPECIMEN`` inside a conda environment with a Python version of 3.10 or higher.

