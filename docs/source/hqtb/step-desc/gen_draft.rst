Step 2: Generate a Draft Model
==============================

Based on the results of step 1, step 2 of the pipelines generates a draft model based on the 
template model.

.. note:: 

    For the pipeline to work, the template model needs to be constructed from the genome that was
    used as a reference during the alignment in step 1.

A graphical representation of the steps can be found below:

.. image:: ../../images/modules/2_gendraft.png

The draft model is constructed based on the idea of :footcite:t:`norsigian2020workflow` by:

- Filtering the BLAST best bidirectional hits by their percentage identity value (PID > threshold suggests, a homolog was found).
- Choosing a medium for the new model. Options include:

    - Using the one from the template
    - Setting all exchanges to open (not advised)
    - Loading one from the ``refineGEMs`` database

- Removing genes from the template that have no (found) homolog in the template genome.

    - Under the condition, that the removed gene is not neccessary for the growth of the model

- Renaming the homologous genes to match the names of the new genome
- Checking unchanged genes

    - remove them, if they are not essential for the growth of the model 

.. footbibliography:: 
