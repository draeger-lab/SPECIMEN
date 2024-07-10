Installation
==============

Installation via github
-----------------------
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

.. hint::

    It is advised to install ``SPECIMEN`` inside a conda environment with a Python version of 3.10 or higher.

