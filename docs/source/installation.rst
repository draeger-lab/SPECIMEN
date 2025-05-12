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

- `DIAMOND, version 2.0.4 or higher <https://github.com/bbuchfink/diamond>`_  (cannot be installed 
  via conda on Windows)

Afterwards, ``SPECIMEN`` can either be accessed via the command line or via importing the package into a Python script.

.. hint::

    The minimum Python requirement for ``SPECIMEN`` is 3.10.

.. warning:: 

    Currently, COBRApy is not yet updated on its main branch for Python 3.13.\\
    Therefore, the highest recommended Python version is 3.12.
    

