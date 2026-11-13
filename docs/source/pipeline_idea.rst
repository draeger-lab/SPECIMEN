:icon:`fa-solid fa-dna` PGAB: From genome sequence to draft model
=================================================================

PGAB stands for Prokaryotic Annotation Pipeline (PGAP)-based workflow.

.. caution::

    This workflow is under development on the ``pgab-dev`` branch and is planned
    for a future update. The overview below describes the planned workflow.

Starting from a genome sequence, PGAB aims to annotate the genome and extend the
protein annotations with database identifiers before generating a draft model. 
After each of the major stages (highlighted in green in :numref:`pgab_workflow`) 
the workflow can be stopped and only the according output will be returned.
:numref:`pgab_workflow` shows the following stages:

1. **Taxonomy check and/or genome annotation:** Uses NCBI's Prokaryotic Genome
   Annotation Pipeline (PGAP) with the genome FASTA file as input. This stage
   produces an annotated FASTA and GFF file and/or taxonomy information in XML/
   TXT format depending on the run settings (``tax-check = 'only'``: Taxonomy 
   information in XML/TXT format, ``tax-check = 'continue'``: All files, 
   ``tax-check = 'none'``: Annotated FASTA and GFF file).
2. Depending on the result from stage 1 and the user-input either route 
   **A - RefSeq Lookup** or route **B - Annotation extension** is selected. If the 
   taxonomy check was performed and PGAP found a taxon with 100% identity and 
   the user did not specify that DIAMOND should be run route A is chosen. In all 
   other cases the next stage is a DIAMOND search to extend the annotation content.

   A. **RefSeq Lookup:** Based on the Genbank identifier with 100% ANI the 
      corresponding RefSeq identifier is searched.
   B. **Annotation extension:** DIAMOND is used to obtain additional protein 
      annotations. The DIAMOND route includes taxonomy-specific and unrestricted 
      searches against NCBI's ``nr`` database, SwissProt, or a user-provided 
      database. The taxonomy-specific search is only performed if a taxid was 
      found with PGAP or if the user provided a taxid.
3. **Draft generation:** Uses either the data from NCBI based on the looked up 
   RefSeq identifier (result from 2.A) or the annotated files from PGAP along 
   with the extended annotations if route 2.B was used with either the 
   :ref:`CMPB<CMPB Overview>` or :ref:`HQTB<HQTB Overview>` workflow to generate 
   the draft model.
4. **Output:** The output varies depending on the following cases:

   A. Only PGAP was requested: The output of PGAP is returned depending on the 
      chosen ``tax-check`` setting.

      .. important::

         **Special case:** ``tax-check = 'only'`` returns only the taxonomy report 
         and never performs any of the other steps.

   B. Route 2.A was used but no draft model generation requested: Only a warning 
      is logged and the results from the steps before are returned.
   C. Route 2.B was used but no draft model generation requested: Only a 
      mapping table in CSV format linking the annotated proteins to database 
      identifiers is returned along with the files from the previous steps.
   D. The user specified that a draft model is generated: The draft model is 
      returned in SBML format and a mapping table in CSV format as well as the 
      files from the previous steps.

.. _pgab_workflow:

.. figure:: images/pgab_pipeline-overview.png
  :alt: Planned PGAB workflow with PGAP annotation, RefSeq lookup and DIAMOND searches, CMPB or HQTB draft generation, and CSV and SBML outputs

  Planned PGAB workflow from a genome sequence to a draft model and a mapping table.
