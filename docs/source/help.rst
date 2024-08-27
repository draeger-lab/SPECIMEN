Help & FAQ
==========

If you run into any problems regarding ``SPECIMEN``
or have wishes and ideas for features to be included, please write an issue on our
`GitHub page <https://github.com/draeger-lab/SPECIMEN>`__ or contact us directly.

Known Bugs and Issues
-----------------------

Pydantic
^^^^^^^^

Pydantic warning :code:`underscore_attrs_are_private has been removed` has not - yet - caused any issues.
However, the core of the problems (= what causes the warning) has yet to be identified. 

BOFdat
^^^^^^

For `BOFdat <https://github.com/jclachance/BOFdat>`_, to be running as expected, 
one needs to change :code:`solution.f` to :code:`solution.objective_value` in the :code:`coenzymes_and_ions.py` file of the tool's files.

.. warning:: 

    The change in the files for BOFdat is neccessary, otherwise it will raise an error while running the connected functions.


FAQs
----

Which types of organism work with ``SPECIMEN``?
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

``SPECIMEN`` was originally written for prokaryota. However, adapting ``SPECIMEN`` to work with 
other organism types is something we hope to archieve with a future update. 
Currently, the workflows have yet to be tested on types of organism other than prokaryota.

Which namespaces can I use?
^^^^^^^^^^^^^^^^^^^^^^^^^^^

While we are working towards namespace independency, currently,
it is advised to stick to the BiGG namespace.
