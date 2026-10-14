Step 2: Generate a Draft Model
==============================

Based on the results of step 1 and the template model, step 2 of the workflow generates a draft model.

.. note:: 

    When running this step outside the context of the workflow, the input files 
    need to be related to each other. The identifiers of the genes 
    for the matches of the bidirectional blast results and the genes in the template model
    need to be the same or can be adjusted using the given parameters as needed.

A graphical representation of the steps can be found below:

.. image:: ../../images/modules/2_gendraft.png

The draft model is constructed based on the idea of :footcite:t:`norsigian2020workflow` by:

1. Filtering the BLAST best bidirectional hits by their percentage identity value (only matches with a PID greater than a given threshold are considered homologous).
2. Choosing a medium for the new model. Options include:

    - Using the one from the template
    - Setting all exchanges to open (not advised, quickly leads to unrealistic growth rates)
    - Loading one from the ``refineGEMs`` database

3. Removing genes from the template that have no (found) homologue in the template genome.

    .. important:: 
        Under the condition, that the removed gene is not necessary for the growth of the model.

4. Renaming the homologous genes to match the names of the new genome
5. Checking unchanged genes

    .. important:: 
        Unchanged genes are removed, if they are not essential for the growth of 
        the model.

        Additionally, genes from the template model, which are part of an enzyme 
        complex, where at least one gene was mapped successfully, are kept as well. 
        They are, however, reported to provide extensive information for possible 
        manual curation after running the workflow.

.. footbibliography:: 
