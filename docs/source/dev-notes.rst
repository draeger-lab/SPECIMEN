Notes for Developers
=====================

To maintain or extend the toolbox ``SPECIMEN`` please install the package via GitHub.

.. hint::

   For help and information about known bugs, refer to :ref:`Help & FAQ`.

Current Status
--------------
Below you can find information on the currently available branches.

If you plan to add a new feature or extension, please create your own branch, as to not implicate other
developers.

``main``
^^^^^^^^

Our first development release of ``SPECIMEN``. Please feel free to test and leave issues on GitHub if you encounter bugs or have suggestions.
The version in this release has as of yet only be tested on prokaryote genomes.

.. note::

    The workflows have yet to be tested on more diverse species and specifically Eukarya, if they work without problems on those as well.
    (Due to the current limitations of the compartments, Eukarya will most likely not work in most cases.)

``docs-update``: ONGOING
^^^^^^^^^^^^^^^^^^^^^^^^

| Documentation branch.
| Place to update, extend, etc. the documentation.
| DO NOT use this branch to edit other folders except for docs! 


``dev``: ONGOING
^^^^^^^^^^^^^^^^

| Development branch. 
| Place to fix errors, chase bugs and polish the code for new releases.
| Additionally, features marked as future update throughout the documentation are actively being worked on.
| If a bug/error needs to be fixed on main, please create a hotfix-branch.


Further branches
^^^^^^^^^^^^^^^^

.. Currently, there are no more branches.

- ``pgab-dev``: Branch for the development of the PGAB workflow 

  - Status: Ongoing work (Master's thesis)


Installation for developers
---------------------------

Into a Conda environment
^^^^^^^^^^^^^^^^^^^^^^^^

Setup a conda virtual environment and use its pip to install ``refineGEMs`` into that environment.
To run the code below, replace ``<EnvName>`` with the name of your environment and ``<Specific Python version >= 3.10>`` with the specific version of Python you want to use.

.. code:: console

   # clone or pull the latest source code
   git clone https://github.com/draeger-lab/specimen.git
   cd specimen

   conda create -n <EnvName> python=<Specific Python version >= 3.10>

   conda activate <EnvName>

   # check that pip comes from <EnvName>
   which pip

   pip install .

This will install all packages denoted in `pyproject.toml`. 

If `which pip` does not show pip in the conda environment you can also create a local environment for which you can 
control the path and use its pip:

.. code:: console

   conda create --prefix ./<EnvName>

   conda activate <path to EnvName>

   <EnvName>/bin/pip install .

Into a Pipenv environment
^^^^^^^^^^^^^^^^^^^^^^^^^

You can use `pipenv <https://pipenv.pypa.io/en/latest/>`__ to keep all dependencies together. Therefore, you will need 
to install ``pipenv`` first. To install ``specimen`` locally complete the following steps:

.. code:: console

   # install pipenv using pip
   pip install pipenv

   # clone or pull the latest source code
   git clone https://github.com/draeger-lab/specimen.git
   cd specimen

   # install local package and dependencies into a virtual environment
   pipenv install .

   # initiate a session in the virtual environment
   pipenv shell

The ``pipenv`` package can also be installed via Anaconda (recommended
if you are a Windows user).

.. hint::

  If you want to be able to safely import the package from anywhere while also retaining the possibility to edit the 
   code, it is recommended to change the :code:`pip install` line from the code blocks to 
   :code:`pip install -e . --config-settings editable_mode=strict`.

Additional packages required for development
--------------------------------------------

.. attention::
    The following packages need to be installed to be able to add content to the `SPECIMEN` documentation.
    
    * `accessible-pygments`
    * `ipython`
    * `nbsphinx`
    * `pandoc`
    * `sphinx`
    * `sphinx_copybutton`
    * `sphinx_rtd_theme`
    * `sphinxcontrib-bibtex`
    * `sphinx-icon`
    

    In addition, `pip-compile` should be installed to update the `requirements.txt` for the next release.

Installing the packages
^^^^^^^^^^^^^^^^^^^^^^^
You can install the packages via pip to your local environment:

.. code:: console
    :class: copyable

    pip install accessible-pygments ipython nbsphinx pandoc sphinx sphinx_copybutton sphinx_rtd_theme sphinxcontrib-bibtex sphinx-icon

.. code:: console
    :class: copyable

    python -m pip install pip-tools

Alternatively, install the tool with the extra `docs`, e.g. 

.. code:: console
    :class: copyable

     pip install -e ".[docs]" --config-settings editable_mode=strict

Todo Tree extension
-------------------

If you are working with VS Code or similar, you can install the ``Todo Tree`` extension and copy the 
content of the ``TodoTree_params.txt`` file in the ``dev`` folder into the corresponding place in your setting 
to enable highlighting and tracing of the ``@KEYWORD`` labels for bugs, discussions and more.

The usage of these keywords strongly recommended, as it make communication between the developers and tracing of 
issues and ideas much easier. The following labels are currently supported:

.. list-table:: Label supported with the SPECIMEN Todo Tree file
    :header-rows: 1
    :widths: 25 75 

    * - label 
      - usage
    * - @TODO
      - something needs to be implemented or changed
    * - @BUG 
      - something is not working as expected, the cause has yet to be determined
    * - @FIXME 
      - known error or issue that needs attention
    * - @TEST 
      - the following code needs to be tested 
    * - @DEPRECATE 
      - the following code can be removed in the future
    * - @NOTE 
      - notes from one dev to another or just a reminder 
    * - @DISCUSSION 
      - an issue, idea or feature, that requires discussion
    * - @DEBUG 
      - label for a debugging switch or similar, see next section
    * - @ASK 
      - when something requires research or further input, before it can be discussed
    * - @IDEA
      - write down ideas for new features, better implementations and more
    * - @WARNING
      - if something needs to be kept in mind or can easily break down, without being a bug

.. note:: 

    For the label to be recognised correctly, the following format is required: ``# @KEYWORD``.


Developer Information
---------------------

More developer-relevant information can be found in the :code:`dev` folder on the `github page <https://github.com/draeger-lab/SPECIMEN>`__.

Updating the `requirements.txt`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
| To create the `requirements.txt` adjust the `requirements.in` file as needed in the folder docs.
| Then navigate to the folder docs in the command line:

.. code-block:: console
    :class: copyable

    cd docs

and use the following command to automatically generate the new `requirements.txt`:

.. code-block:: console
    :class: copyable
    
    python3 -m piptools compile --strip-extras --output-file=requirements.txt requirements.in

To bump to the newest versions possible, use the following command in the `docs` directory:

.. code-block:: console
    :class: copyable

    pip-compile --upgrade


Debugging switches
^^^^^^^^^^^^^^^^^^

- You can enable debug logging by replacing ``level=logging.INFO``  with ``level=logging.DEBUG``.
- If you want your print message to show in the log file, replace the ```print()`` statement by ``logging.info()``.
- For debugging of pandas warnings or issues ``pd.options.mode.chained_assignment = None`` needs to be commented out.
- | Additionally, some modules contain comment blocks inf the format shown below.
  | By enabling the code lines between the dotted lines, a debugging-mode is run, which e.g. subsets the data to shorten the runtime to make debugging faster.

.. code-block:: python
    :class: copyable

    # @DEBUG ...............
    # some code
    # ......................


Documentation Notes
-------------------

The documentation is generated based on the Sphinx :code:`sphinx.ext.autodoc` extension.
A mustache-file with additional formatting can be found in the :code:`dev` folder (ready to integrate in e.g. VSCode). 
For further information refer to the `refineGEMS documentation notes <https://refinegems.readthedocs.io/en/latest/development.html>`__.

Please annotate new functions with restructured text for them to be included into this documentation.
Furthermore, please use type hinting for the functions, e.g.:

.. code-block:: python

    def func(num:int, square:bool=True) -> int:
        ...
        return x