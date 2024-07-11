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

    The pipeline has yet to be tested on different species and specifically Eukarya, if it works without problems on those as well.

``dev``: ONGOING
^^^^^^^^^^^^^^^^

Development branch. 
Currently, features marked as future update throughout the documentation are worked on.

Further branches
^^^^^^^^^^^^^^^^

Currently, there are no more branches.

Developer Information
---------------------

More developer-relevant information can be found in the :code:`dev` folder on the `github page <https://github.com/draeger-lab/SPECIMEN>`__.

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