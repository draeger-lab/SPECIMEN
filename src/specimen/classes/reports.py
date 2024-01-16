"""Classes to generate, handle, manipulate and save reports.
"""

__author__ = 'Carolin Brune'

################################################################################
# requirements
################################################################################

import cobra
from importlib.resources import files
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import pandas as pd
import seaborn as sns
import warnings

from pathlib import Path

from . import medium

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


class GrowthAnalysisReport():
    """GrowthAnalysisReport saved the results of analyse_growth() and allows for saving and visualising the results.

    :param growth_sim_results: The results of  :py:func:`specimen.classes.medium.growth_simulation`
    :type growth_sim_results: dict, medium name : growth rate (optional)
    :param minimal_medium: The minimal medium for the underlying model.
    :type minimal_medium: medium.Medium (optional)
    :param auxo_results: The results of auxotrophies_simulation().
    :type auxo_results: nested dict

    :var growth_rates: The results of  medium.growth_simulation()
    :vartype growth_rates: dict, medium name : growth rate (optional)
    :var minimal_medium: The minimal medium for the underlying model.
    :vartype minimal_medium: medium.Medium (optional)
    :var auxotrophy_results: The results of auxotrophies_simulation().
    :vartype auxotrophy_results: nested dict
    """

    def __init__(self, growth_sim_results=None, minimal_medium=None, auxo_results=None):

        self.growth_rates = growth_sim_results
        self.minimal_medium = minimal_medium
        self.auxotrophy_results = auxo_results


    def growth_rates_table(self):
        """Turn the dictionary of growth rates into a DataFrame.

        :returns: The DataFrame.
        :rtype: pd.DataFrame or None
        """

        if self.growth_rates:
            return pd.DataFrame({'medium':self.growth_rates.keys(),
                            'growth':self.growth_rates.values()})
        else:
            return


    def visualise_growth(self, unit='h'):
        """Visualise thegrowth of the simulation results.
        Either the growth rates themselves can be visualised per hour
        or per min or the doubling times in minutes.

        :param unit: Determines the type of visualisation that is prodoced.
            Can be 'h', 'min' or 'dt'. First two ones are the unit for the
            growth rate visualisation, the last one makes the doubling time graphic.
            Default is 'h'.
        :type unit: string, choices ['h','min','dt']

        :returns: The visualisation.
        :rtype: plt.fig
        """

        data = self.growth_rates_table()
        # create y-axis data and xlabels for the given unit
        match unit:
            case 'h':
                ydata = data['growth']
                xlab = r'growth rate $[\frac{mmol}{gDWh}]$'
            case 'min':
                ydata = data['growth']/60
                xlab = r'growth rate $[\frac{mmol}{gDWmin}]$'
            case 'dt':
                ydata = 70/((data['growth']/60)*100)
                xlab = 'doubling time [min]'
            case _:
                raise ValueError(F'Unknown unit for visualise growth: {unit}')

        # re-check data values
        ydata = [x if x>0.0005 else 0.0 for x in ydata]

        # create a colour gradient
        cmap = sns.cubehelix_palette(start=0.5, rot=-.75, gamma=0.75, dark=0.25, light=0.6, reverse=True, as_cmap=True)
        # set up the figure
        fig = plt.figure()
        ax = fig.add_axes([0,0,1,1])
        # construct the plot
        max_ydata = max(ydata)
        if max_ydata <= 0:
            warnings.warn('Model is not able to grow on all media. Returning empty figure.')
            return fig
        cont = ax.barh(data['medium'], ydata, color=cmap([_/max_ydata for _ in ydata]),
               label=data['medium'])
        ax.bar_label(cont, fmt='%.2f', color='black', padding=1.0)
        ax.set_ylabel('medium', labelpad=12)
        ax.set_xlabel(xlab, labelpad=12)
        xlims = ax.get_xlim()
        ax.set_xlim(xlims[0],xlims[1]+(xlims[1]*0.03))
        ax.set_title(F'Simulated growth results')

        return fig


    def auxotrophies_table(self):
        """Tranform the auxotrophy_results attribute into a table format.

        :returns: the table (amino acids vs medium names)
        :rytpe: pd.DataFrame
        """
        if self.auxotrophy_results:
            table = pd.DataFrame.from_dict(self.auxotrophy_results)
            return table
        else:
            return pd.DataFrame()


    def visualise_auxotrophies(self):
        """Visualise the results of the auxotrophy tests as a heatmap.

        :returns: The generate heatmap.
        :rtype: plt.Figure or None
        """

        if self.auxotrophy_results:

            # create table
            data = self.auxotrophies_table()
            # create a colour gradient
            cmap = sns.cubehelix_palette(start=0.5, rot=-.75, gamma=0.75, dark=0.25, light=0.6, reverse=True, as_cmap=True)
            # set up the figure
            fig = plt.figure()
            ax = fig.add_axes([0,0,1,1])
            # # create heatmap
            sns.heatmap(data, ax=ax, cmap=cmap, cbar_kws={'label': 'flux'}, annot = True, fmt='.2f')
            # add labels
            ax.set_ylabel('amino acid', labelpad=12)
            ax.set_xlabel('medium', labelpad=12)
            ax.set_title('Fluxes for auxotrophy tests')

            return fig
        else:
            return


    def save(self,dir):
        """Save the report.

        :param dir: Path to a directory to save the output to.
        :type dir: string
        """

        # make sure given directory path ends with '/'
        if not dir.endswith('/'):
            dir = dir + '/'

        # collect all produced file in one directory
        try:
            Path(F"{dir}growth-analysis/").mkdir(parents=True, exist_ok=False)
            print(F'Creating new directory {F"{dir}growth-analysis/"}')
        except FileExistsError:
            print('Given directory already has required structure.')


        if self.growth_rates:
            # save the visualisation of the growth rates
            growth_vis = self.visualise_growth('h')
            growth_vis.savefig(F'{dir}growth-analysis/visualise_reactions_h.png', bbox_inches='tight')
            growth_vis = self.visualise_growth('min')
            growth_vis.savefig(F'{dir}growth-analysis/visualise_reactions_min.png', bbox_inches='tight')
            # visualise doubling time
            growth_vis = self.visualise_growth('dt')
            growth_vis.savefig(F'{dir}growth-analysis/visualise_reactions_dt.png', bbox_inches='tight')

            # save the growth rates as tabular information
            self.growth_rates_table().to_csv(F'{dir}growth-analysis/table_growth_rates.tsv', sep='\t', index=False)

        if self.minimal_medium:
            medium.save_db({self.minimal_medium.name:self.minimal_medium}, F'{dir}growth-analysis/minimal_medium.csv')

        if self.auxotrophy_results:
            # save the visualisation of the growth rates
            aux_vis = self.visualise_auxotrophies()
            aux_vis.savefig(F'{dir}growth-analysis/visualise_auxotrophies.png', bbox_inches='tight')

            # save the growth rates as tabular information
            self.auxotrophies_table().to_csv(F'{dir}growth-analysis/auxotrophies_rates.tsv', sep='\t', index=True)


