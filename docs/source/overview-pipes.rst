Overview of the Workflows
=========================

The following workflows are currently available:

- ``CMPB``: CarveMe + ModelPolisher based 
- ``HQTB``: High-quality template based 

More information about these differnt types of workflows can be found below.

.. hint:: 

    Which type of workflow to use is heavily dependent on your 
    organism, available information for your organism 
    and the data you wish to use as input.

CarveMe + ModelPolisher based (``CMPB``) workflow
-------------------------------------------------

The ``CMPB`` workflow is based on generating a draft model using `CarveMe <https://github.com/cdanielmachado/carveme>`__ and subsequently extending and polishing
the model using various tools.

This workflow can be run on an annotated genome or an already generated CarveMe model and requires very little additional
information to be run on its base settings. 
However, additional information can be added from 
e.g. KEGG or BioCyc to perform an automated gap filling using `refineGEMs <https://github.com/draeger-lab/refinegems/tree/dev-2>`__ :footcite:p:`bauerle2023genome`.

.. note::

    Currently, the ModelPolisher connection is still under construction, but the workflow
    can already be run.

.. toctree::
    :maxdepth: 2
    :caption: Further Information 

    About CMPB <cmpb/about-pipeline>
    Run CMPB <cmpb/run-pipeline>
    CMPB Configuration <cmpb/cmpb-config>

High-quality template based (``HQTB``) workflow
-----------------------------------------------

.. warning:: 
    Due to chances in ``refineGEMs``, this workflow is under heavy 
    developement and may not work as expected.

The ``HQTB`` workflow curates a new model from an annotated genome based on a high-quality template model 
(plus corresponding annotated genome) and additional database information. 

This workflow aims to profit from already performed (manual) curation of an already existing model, 
to carry this knowledge into the new model. The closer the template is to the original, the more knowledge 
can potentially be carried over. Therefore, this workflow is more useful, if the user already has a model of
a similar organism compared to the one for which the new model should be curated for.

.. toctree::
    :maxdepth: 2
    :caption: Further Information 

    About HQTB <hqtb/about-pipeline>
    Run HQTB <hqtb/run-pipeline>
    HQTB Configuration <hqtb/hqtb-config>

More ideas for workflows
------------------------

Below are some ideas for workflows, to be implemented in future update(s):

.. toctree::
    PGAB <pipeline_idea>