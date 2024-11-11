Notes for Developers
=====================

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

    The workflow has yet to be tested on different species and specifically Eukarya, if it works without problems on those as well.

``dev``: ONGOING
^^^^^^^^^^^^^^^^

| Development branch. 
| Place to fix errors, chase bugs and polish the code.
| Additionally, features marked as future update throughout the documentation are actively being worked on.

Further branches
^^^^^^^^^^^^^^^^

Currently, there are no more branches.

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