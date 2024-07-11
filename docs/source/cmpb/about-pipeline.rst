About ``CMPB``
==============

The CarveMe + ModelPolisher based (CMPB) pipeline curates a model using refineGEMs and ModelPolisher.
The starting point is either the input files for CarveMe (future update) or an already built model. If the model is/was 
built with CarveMe a correction is performed. After that the model is gap filled. The gap fill step can be either done 
with a manually created file or automatically. Afterwards ModelPolisher is used to enhance the annotation content 
(future update). Further annotation enhancement is done by adding Pathways as Groups from KEGG and using the SBOannotator 
via the refineGEMs connection to get more specific SBO term annotations. Afterwards the model is cleaned up with 
MassChargeCuration and, optionally, with BOFdat. Both tools are accessed via refineGEMs. A further option for the clean 
up is the removal of duplicates. There are optional steps to analyse the model with FROG (future update) or MEMOTE. The 
default analysis steps include model statistics, growth analysis and auxotrophy tests. For each step the model version and, 
optionally, the according MEMOTE report can be saved.

Overview of the CMPB Pipeline
-----------------------------

The following image shows an overrview of the steps of the pipeline:

.. _cmpb_workflow:

.. figure:: ../images/cmpb_pipeline-overview.png
  :alt: Workflow from CaveMe to close-to-final model

  Workflow from CaveMe to close-to-final model
