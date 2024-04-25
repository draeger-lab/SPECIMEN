Notes for Developers
=====================

Current Status
--------------
Below you can find information on the currently available branches.

If you plan to add a new feature or extension, please create your own branch, as to not implicate other
developers.

**Main** 

The current version of the pipeline has been successfully tested on a *Klebsiella pneumoniae* and *Klebsiella oxytoca* genomes using
a *Klebsiella sp.* reference database and a *Klebsiella pneumoniae* template model.

.. note::

    The pipeline has yet to be tested on different species and specifically Eukarya, if it works without problem on those as well.

**refinegems-integration** : ONGOING

Branch for the integration of SPECIMEN-functionalities into refineGEMS and vice versa.

Integration almost done. Cleanup for first release of SPECIMEN pending

**Further branches**

Currently, there are no more branches.

Developer Information
---------------------

More developer-relevant information can be found in the :code:`dev` folder on the `github page <https://github.com/draeger-lab/SPECIMEN>`__.

Documentation Notes
-------------------
The documentation is generated based the Sphinx :code:`sphinx.ext.autodoc` extension.
A mustache-file with additional formatting can be found in the :code:`dev` folder (ready to integrate in e.g. VSCode). 
For further information refer to the `refineGEMS documentation notes <https://refinegems.readthedocs.io/en/latest/development.html>`__.

Please annotate new functions with restructured text for them to be included into this documentation.
Furthermore, please use type setting for the functions, e.g.:

.. code-block:: python

    def func(num:int, square:bool=True) -> int:
        ...
        return x