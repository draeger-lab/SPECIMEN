Help & FAQ
==========

If you run into any problems regarding **SPECIMEN** 
or have wishes and ideas for features to be included, please write an issue on our
`github page <https://github.com/draeger-lab/SPECIMEN>`__ or contact us directly.

Known Bugs and Issues
-----------------------

Pydantic
^^^^^^^^

Pydantic warning `underscore_attrs_are_private has been removed` has not - yet - caused any issues.
However, the core of the problems (= what causes the warning) has yet to be identified. 


FAQs
----

Which types of organism work with ``SPECIMEN``?
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

``SPECIMEN`` was originally written for prokaryota. However, adapting ``SPECIMEN`` to work with 
other organism types is something we hope to archieve with a future update. 
Currently, the pipelines have yet to be tested on types of organism other than prokaryota.

Which namespaces can I use?
^^^^^^^^^^^^^^^^^^^^^^^^^^^

While we are working towards namespace independency, currently,
it is advised to stick to the BiGG namespace.
