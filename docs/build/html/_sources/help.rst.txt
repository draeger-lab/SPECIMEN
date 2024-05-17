Help & FAQ
==========

If you run into any problems reagrding **SPECIMEN** 
or have wishes and ideas for features to be included, please write an issue on our
`github page <https://github.com/draeger-lab/SPECIMEN>`__ or contact us directly.

Known Bugs and Issues
-----------------------

Pydantic
^^^^^^^^

Pydantic warning `underscore_attrs_are_private has been removed` has not - yet - caused any issues
however, the core of the problems (= what causes the warning) has yet to be identifies. 


FAQs
----

Which organisms work with ``SPECIMEN``?
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

``SPECIMEN`` was originally written for prokaryota, however adapting to work with 
other organisms is something we hope to archieve witha future update. 
Currently, the pipelines have yet to be tested on organisms other than prokaryota. 

Which namespaces can I use?
^^^^^^^^^^^^^^^^^^^^^^^^^^^

While we are working towards namespace indepancy of this pipeline, currently,
it is advised to stick to the BiGG namespace.