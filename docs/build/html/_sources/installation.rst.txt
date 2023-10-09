Installation
==============

Installation via github
-----------------------
Download the **SPECIMEN** repository `here <https://github.com/cb-Hades/klebsiella-pipeline/tree/main>`_ and run the command :code:`pip install -e .` inside the top-level directory.

For all options to run error-free, the following tools need to be installed additionally:

- `MCC - MassChargeCuration <https://github.com/Biomathsys/MassChargeCuration/tree/main/MCC>`_
- `DIAMOND, version 2.0.4 or higher <https://github.com/bbuchfink/diamond>`_
- `EntrezDirect <https://www.ncbi.nlm.nih.gov/books/NBK179288/>`_, if no NCBI mapping has been created beforehand
- `BOFdat <https://github.com/jclachance/BOFdat>`_, for running as expected, one needs to change :code:`solution.f` to :code:`solution.objective_value` in the :code:`coenzymes_and_ions.py` file of the tool's files

Afterwards, the **SPECIMEN** can either be accessed via the command line or via importing the package into a Python script.

.. note::
    It is avised to install **SPECIMEN** inside a conda environment with a python version of 3.10 or higher.
