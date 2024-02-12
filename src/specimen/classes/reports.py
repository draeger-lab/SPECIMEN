"""Classes to generate, handle, manipulate and save reports.
"""

__author__ = 'Carolin Brune'

################################################################################
# requirements
################################################################################

import pandas as pd

from refinegems.classes.reports import ModelInfoReport

################################################################################
# variables
################################################################################

################################################################################
# classes
################################################################################

class SpecimenModelInfoReport(ModelInfoReport):
    
    def __init__(self, model) -> None:

        # call the superclass
        super().__init__(model)

        # find out the origin of the reactions
        reac_origin_counts = {'via template':0, 'via MetaNetX':0, 'via KEGG':0, 'via gapfilling':0, 'else':0}
        for reac in model.reactions:
            # get origin of reaction (based on workflow notation)
            if 'creation' in reac.notes.keys():
                if reac.notes['creation'] in reac_origin_counts.keys():
                    reac_origin_counts[reac.notes['creation']] += 1
                else:
                    reac_origin_counts['else'] += 1
            else:
                reac_origin_counts['else'] += 1
        
        # add new attribute
        self.reac_origin_c = reac_origin_counts

    # extemd format table function from parent class
    def format_table(self) -> pd.DataFrame:
        table = super().format_table()
        table['#reaction origin'] = str(self.reac_origin_c).replace('{',r'').replace('}',r'').replace('\'',r'')
        return table
    
    # depending on the implementation, save and make html 
    # can be inherited or need to be overwritten 

