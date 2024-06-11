Overview of the Pipelines
=========================

The following pipelines are currently available:

- ``CMPB``: CarveMe + ModelPolisher based 
- ``HQTB``: High-quality template based 

More information about these differnt types of pipelines can be found below.

.. hint:: 

    Which type of pipeline to use with out data is heavily dependant on your organism, 
    available information for it and the data you wish to use as input.

CarveMe + ModelPolisher based (CMPB) pipeline
---------------------------------------------

The CMPB pipeline is based on generating a draft model using CarveMe and subsequently extending and polishing
the model using various tools.

This pipeline can be run on an annotated genome or an already generated CarveMe model and requires very little additional
information to be run on its base settings. 
However additional information can be added from 
e.g. KEGG or BioCyc to perform an automated gapfilling using `refineGEMs <https://github.com/draeger-lab/refinegems>`__ :footcite:p:`bauerle2023genome`.

.. note::

    Currently, the ModelPolisher connections is still under construction, but the pipeline
    can already be run.

.. toctree::
    :maxdepth: 2
    :caption: Further Information 

    About CMPB <cmpb/about-pipeline>
    Run CMPB <cmpb/run-pipeline>
    CMPB Configuration <cmpb/cmpb-config>

High-quality template based (HQTB) pipeline
-------------------------------------------

The HQTB pipeline curates a new model from an annotated genome based on a high-quality template model 
(plus corresponding annotated genome) and additional database information. 

This pipeline aims to profit from already performed (manual) curation of the already existing model, 
to carry this knowledge into the new model. The closer the template is to the original, the more knowledge 
can potential be carried over. Therefore, this pipeline is more useful, if the user already has a model of
a similar organism to the one the new model should be curated for.

.. toctree::
    :maxdepth: 2
    :caption: Further Information 

    About HQTB <hqtb/about-pipeline>
    Run HQTB <hqtb/run-pipeline>
    HQTB Configuration <hqtb/hqtb-config>

More ideas for pipelines
------------------------

Below are some ideas for pipelines, which may be implemented in future update(s):

.. toctree::
    PGAB <pipeline_idea>