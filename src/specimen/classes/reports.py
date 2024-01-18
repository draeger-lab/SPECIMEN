"""Classes to generate, handle, manipulate and save reports.
"""

__author__ = 'Carolin Brune'

################################################################################
# requirements
################################################################################

from pathlib import Path

################################################################################
# variables
################################################################################

################################################################################
# classes
################################################################################

class ModelStatisticsReport():
    """ModelStatisticsReport objects are used to store statistical information
    about a model. Additionally, this class provides functionalities to report
    said information.

    :param model_id: The model's ID (its name).
    :type model_id: string
    :param total_reac: Total number of reactions in the model.
    :type total_reac: int
    :param total_gene: Total number of genes in the model.
    :type total_gene: int
    :param total_meta: Total number of metabolites in the model.
    :type total_meta: int
    :param reac_origin_counts: Descriptors and corresponding count for where the reactions where build from.
        Should be a dictionary of string and ints.
    :type reac_origin_counts: dict
    :param reac_with_gpr: Number of reactions that are associated with a gene reaction rule (GPR)
    :type reac_with_gpr: int

    :var model_id: The model's ID the statistic is about.
    :vartype model_id: string
    :var reac: Information about the reaction objects, including "total" count, "origin" count and number of reactions "with gpr".
    :vartype reac: dict (string as above, int or dict as value)
    :var meta: Information about the metabolite objects, including "total" count.
    :vartype meta: dict (string as above, int as value)
    :var gene: Information about the gene objects, including "total" count.
    :vartype gene: dict (string as above, int as value)
    """

    def __init__(self, model_id, total_reac, total_meta, total_gene,
                 reac_origin_counts, reac_with_gpr):

        self.model_id = model_id
        self.reac = {'total': total_reac, 'origin':reac_origin_counts, 'with gpr':reac_with_gpr}
        self.meta = {'total': total_meta}
        self.gene = {'total': total_gene}

    def get_statistics(self):
        """Generate a string got the statistics.

        @TODO find a better format, maybe html?

        :returns: the statistics
        :rtype: string
        """

        report = F'''Statistical report for model {self.model_id}

reactions
---------
total : {self.reac['total']}
origin:
    {(chr(10)+"    ").join([F"{k}: {v}" for k,v in self.reac['origin'].items()])}
with gpr : {self.reac['with gpr']}

metabolites
-----------
total : {self.meta['total']}

genes
-----
total : {self.gene['total']}
'''
        return report


    def save(self, dir):
        """Save the report.

        :param dir: Path to a directory to save the output to.
        :type dir: string
        """
        # make sure given directory path ends with '/'
        if not dir.endswith('/'):
            dir = dir + '/'

        # collect all produced file in one directory
        try:
            Path(F"{dir}statistics/").mkdir(parents=True, exist_ok=False)
            print(F'Creating new directory {F"{dir}statistics/"}')
        except FileExistsError:
            print('Given directory already has required structure.')

        # save the statistics report
        with open(F'{dir}/statistics/statistics_report.txt','w') as f:
            f.write(self.get_statistics())

