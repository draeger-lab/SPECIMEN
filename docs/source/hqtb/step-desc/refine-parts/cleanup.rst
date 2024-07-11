Step 3, Part 2: Cleaning the Model
==================================

The third part of the refinement aims to clean the model regarding its entities, including:

- Checking the directionality of reactions using BioCyc.
- Complete BioCyc/MetaCyc, if only one of the two is present.
- Find and resolve duplicate metabolites and reactions.
- Delete unused metabolites.
- Perform (additional) gap filling.

Except for the *complete Bio/MetaCyc annotations* all steps are optional.

.. image:: ../../../images/modules/3_2_cleanup.png

.. warning::

    The gap filling is currently only available in the COBRApy variant.

    This part of the pipeline is still a working process, 
    stay tuned for future updates.

