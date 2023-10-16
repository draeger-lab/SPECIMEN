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

KEGG_GLOBAL_PATHWAY = {'01100': 'Metabolic pathways',
                       '01110': 'Biosynthesis of secondary metabolites',
                       '01120': 'Microbial metabolism in diverse environments'}

KEGG_OVERVIEW_PATHWAY = {'01200': 'Carbon metabolism',
                         '01210': '2-Oxocarboxylic acid metabolism',
                         '01212': 'Fatty acid metabolism',
                         '01230': 'Biosynthesis of amino acids',
                         '01232': 'Nucleotide metabolism',
                         '01250': 'Biosynthesis of nucleotide sugars',
                         '01240': 'Biosynthesis of cofactors',
                         '01220': 'Degradation of aromatic compounds'}

KEGG_METABOLISM_PATHWAY = files('KEGG').joinpath('KEGG_pathway_metabolism.csv')

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


class PanCoreAnalysisReport():
    """PanCoreAnalysisReport objects are used to store, handle and present the
    information from the pan-core analysis (comparison of a model and a pan-core model).

    :param model: The model the test was performed on.
    :type  model: cobra.Model
    :param core_reac: List of reactions characterised as core.
    :type core_reac: list, elements are cobra.Reaction
    :param pan_reac: List of reactions characterised as pan.
    :type pan_reac: list, elements are cobra.Reaction
    :param novel_reac: List of reactions characterised as core.
    :type novel_reac: list, elements are cobra.Reaction

    :var model: The model the test was performed on.
    :vartype  model: cobra.Model
    :var core_reac: List of reactions characterised as core.
    :vartype core_reac: list, elements are cobra.Reaction
    :var pan_reac: List of reactions characterised as pan.
    :vartype pan_reac: list, elements are cobra.Reaction
    :var novel_reac: List of reactions characterised as core.
    :vartype novel_reac: list, elements are cobra.Reaction
    """

    def __init__(self, model,
                 core_reac=None, pan_reac=None, novel_reac=None):

        # general attributes
        self.model = model
        # reaction attributes
        self.core_reac = core_reac
        self.pan_reac = pan_reac
        self.novel_reac = novel_reac
        # ...

    def get_reac_counts(self):
        """Return a dictionary of the counts of the reactions types (core, pan, novel).
        """

        counts = {}
        if self.core_reac:
            counts['core'] = len(self.core_reac)
        else:
            counts['core'] = 0
        if self.pan_reac:
            counts['pan'] = len(self.pan_reac)
        else:
            counts['pan'] = 0
        if self.novel_reac:
            counts['novel'] = len(self.novel_reac)
        else:
            counts['novel'] = 0

        return counts


    #@TODO
    def isValid(self,check='reaction-count'):
        """Check if a certaing part of the analysis is valid.

        Currently possible checks:
            reaction-count : check if the number of reactions in the model
                             equal the sum of the novel, pan and core reactions

        @TODO
            implements more checks

        :param check: Part of the report to check if it is valid.
            Default is 'reaction-count'. Options are in the function description.
        :type check: string

        :raises: :class:`ValueError`: 'Unknown string for parameter check: '

        :returns: True, if test was successful.
        :rtype: bool
        """

        match check:
            case 'reaction-count':
                pc_total = sum(self.get_reac_counts().values())
                diff = len(self.model.reactions) - pc_total
                if diff != 0:
                    return False
                else:
                    return True
            case _:
                raise ValueError('Unknown string for parameter check: ', check)


    def visualise_reactions(self):
        """Visualise the results of the pan-core analysis for the reactions as a donut chart.

        :returns: The generated visualisation object (donut chart).
        :rtype: matplotlib.figure.Figure
        """

        # check the counts
        # ----------------
        counts = self.get_reac_counts()

        if self.isValid('reaction-count'):
            vis_data = list(counts.values())
            vis_label = list(counts.keys())
            vis_color = ['lightgreen','lightskyblue','lightcoral']

        else:
            warnings.warn('Discrepancies between number of reactions in model and sum of novel, pan and core reactions detected.')
            vis_data = list(counts.values()).append(len(self.model.reactions) - sum(counts.values))
            vis_label = list(counts.keys()).append('discrepancies')
            vis_color = ['lightgreen','lightskyblue','lightcoral', 'gold']

        # plot a donut chart
        # ------------------
        fig, ax = plt.subplots()
        wedges, texts, autotexts = ax.pie(vis_data,
                                          autopct=lambda pct: "{:.1f}%\n({:.0f})".format(pct, (pct/100)*sum(vis_data)),
                                          pctdistance=.8, labeldistance=1.25,
                                          colors=vis_color,
                                          radius = 1, wedgeprops=dict(width=0.4, edgecolor='w'))

        ax.legend(wedges, vis_label,
                  title="Classification",
                  loc="center left",
                  bbox_to_anchor=(1, 0, 0.5, 1))

        ax.set_title(F"Results of core-pan analysis\nof the reactions of model {self.model.id}")

        return fig


    #@TODO
    def save(self, dir):
        """Save the results inside a PanCoreAnalysisReport object.

        The function creates a new folder 'pan-core-analysis'
        inside the given directory and creates the following documents:

        - table_reactions.tsv : reactions ID mapped to their labels
        - visualise_reactions : donut chart of the values above

        :param dir: Path to a directory to save the output to.
        :type dir: string
        """
        # ..........................................................
        #@TODO
        #    - an easily human readable overview file?
        #    - smth about metabolites?
        # ..........................................................

        # make sure given directory path ends with '/'
        if not dir.endswith('/'):
            dir = dir + '/'

        # collect all produced file in one directory
        try:
            Path(F"{dir}pan-core-analysis/").mkdir(parents=True, exist_ok=False)
            print(F'Creating new directory {F"{dir}pan-core-analysis/"}')
        except FileExistsError:
            print('Given directory already has required structure.')

        # save the reactions visualisation
        reac_vis = self.visualise_reactions()
        reac_vis.savefig(F'{dir}pan-core-analysis/visualise_reactions.png', dpi=reac_vis.dpi)

        # save table of reactions mapped to characterisation
        if not self.isValid(check='reaction-count'):
            warnings.warn('Discrepancies between number of reactions in model and sum of novel, pan and core reactions detected. Only labbeld reactions will be written to table.')
        reac_tab = pd.DataFrame({'reaction_id': [_.id for _ in self.core_reac] + [_.id for _ in self.pan_reac] + [_.id for _ in self.novel_reac],
                                'pan-core': (['core']*len(self.core_reac)) + (['pan']*len(self.pan_reac)) + (['core']*len(self.novel_reac))})
        reac_tab.to_csv(F'{dir}pan-core-analysis/table_reactions.tsv', sep='\t', index=False)


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


class PathwayAnalysisReport():
    """PathwayAnalysisReport saves the results of kegg_pathway_analysis() (see core.analysis)
    and allows for saving and visualising the results.

    :param total_reac: Total number of reactions in the model
    :type total_reac: int (optional)
    :param kegg_count: Number of reactions with KEGG pathway annotation.
    :type kegg_count: int (optional)
    :param kegg_global: Identifiers and counts for the global level.
    :type kegg_global: dict(str,int)
    :param kegg_over: Identifiers and counts for the overview level.
    :type kegg_over: dict(str,int) (optional)
    :param kegg_rest: Identifiers and counts not match by overview and global level.
    :type kegg_rest: dict(str,int) (optional)

    :var total_reac: Total number of reactions in the model
    :vartype total_reac: int
    :var kegg_count: Number of reactions with KEGG pathway annotation.
    :vartype kegg_count: int
    :var kegg_global: Identifiers and counts for the global level.
    :vartype kegg_global: dict(str,int)
    :var kegg_over: Identifiers and counts for the overview level.
    :vartype kegg_over: dict(str,int)
    :var kegg_paths: Identifiers and counts not match by overview and global level.
    :vartype kegg_paths: dict(str,int)
    """

    def __init__(self,
                 total_reac=None, kegg_count=None,
                 kegg_global=None, kegg_over=None, kegg_rest=None):

        # general counts
        self.total_reac = total_reac
        self.kegg_count = kegg_count

        # kegg pathways
        self.kegg_global = kegg_global
        self.kegg_over = kegg_over
        self.kegg_paths = kegg_rest


    def visualise_kegg_counts(self):
        """Visualise the amounts of reaction with and without
        KEGG pathway annotation.

        :returns: The generated plot.
        :rtype: matplotlib.pyplot.figure
        """

        explode = (0.0,0.1)
        fig, ax = plt.subplots()
        values = [self.kegg_count, self.total_reac-self.kegg_count]
        labels = ['yes','no']
        ax.pie(values,
               autopct=lambda pct: "{:.1f}%\n({:.0f})".format(pct, (pct/100)*sum(values)),
               colors=['lightgreen','lightskyblue'],
               explode=explode, shadow = True, startangle=90)
        ax.legend(labels, title='KEGG\npathway')

        return fig


    def visualise_kegg_pathway(self, plot_type='global', label='id'):
        """Visualise the KEGG pathway identifiers present.

        Depending on the :plot_type:, different levels of pathway identifiers
        are plotted:

        - global: check and plot only the global pathway identifiers
        - overview: check and plot only the overview pathway identifiers
        - high: check and plot all identifiers grouped by their high level pathway identifiers. This option uses label=name, independedly of the input
        - all: check and plot all identifiers

        :param plot_type: Set the type of plot. Can be global, overview, high or existing.
            Default is global.
        :type plot_type: string
        :param id: Set the pathway label to identifier (ID) or actual name.
            Default is 'id', other option is 'name'.
        :type id: string

        :returns: The generated plot.
        :rtype: matplotlib.pyplot.figure
        """

        # get data and KEGG pathway label mapping
        # for the given plot type
        match plot_type:
            case 'global':
                data = self.kegg_global
                label_map = KEGG_GLOBAL_PATHWAY
                title_type = 'global identifiers'
            case 'overview':
                data = self.kegg_over
                label_map = KEGG_OVERVIEW_PATHWAY
                title_type = 'overview identifiers'
            case 'high':
                label = 'name'
                data = self.kegg_paths
                label_map = pd.read_csv(KEGG_METABOLISM_PATHWAY, dtype=str).set_index('id')[['group']].to_dict()['group']
                title_type = 'grouped identifiers'
            case 'existing':
                data = self.kegg_paths
                label_map = pd.read_csv(KEGG_METABOLISM_PATHWAY, dtype=str).set_index('id')[['specific']].to_dict()['specific']
                title_type = 'identifiers in model'
            case _:
                warnings.warn(F'Unknown option for plot_type, choosing "global" istead: {plot_type}')
                data = self.kegg_global
                label_map = KEGG_GLOBAL_PATHWAY
                title_type = 'global identifiers'

        # create the label
        match label:
            case 'id':
                for k in label_map:
                    if k not in data and plot_type != 'existing':
                        data[k] = 0
            case 'name':
                old_data = data
                data = {}
                for k in label_map:
                    if k not in old_data:
                        if plot_type != 'existing' and label_map[k] not in data:
                            data[label_map[k]] = 0
                    else:
                        if label_map[k] not in data:
                            data[label_map[k]] = old_data[k]
                        else:
                            data[label_map[k]] += old_data[k]

            case _:
                warnings.warn(F'Unknown input for label: {label}. Using "id" instead.')
                for k in label_map:
                    if k not in data:
                        data[k] = 0

        data = pd.DataFrame(data.items(), columns=['label', 'counts']).sort_values('counts')
        cdata = data.counts.values
        ldata = data.label.values

        # create the graph
        # create a colour gradient
        cmap = sns.cubehelix_palette(start=0.5, rot=-.75, gamma=0.75, dark=0.25, light=0.6, reverse=True, as_cmap=True)
        # set up the figure
        fig = plt.figure()
        if 'existing' == plot_type:
            ax = fig.add_axes([0,0,5,6])
            # construct the plot
            cont = ax.barh(ldata, cdata, color=cmap([x/max(cdata) for x in cdata]),
                   label=ldata)
            ax.bar_label(cont, fmt='%d', color='black', padding=1.0)
            ax.set_ylabel('KEGG pathway')
            ax.set_xlabel('Number of reaction annotations')
            xlims = ax.get_xlim()
            ax.set_xlim(xlims[0],xlims[1]+(xlims[1]*0.03))
            ax.set_title(F'Pathway analysis with KEGG: {title_type}')
        else:

            ax = fig.add_axes([0,0,1,1])
            # construct the plot
            cont = ax.barh(ldata, cdata, color=cmap([x/max(cdata) for x in cdata]),
                   label=ldata)
            ax.bar_label(cont, fmt='%d', color='black', padding=1.0)
            ax.set_ylabel('KEGG pathway', labelpad=12)
            ax.set_xlabel('Number of reaction annotations', labelpad=12)
            xlims = ax.get_xlim()
            ax.set_xlim(xlims[0],xlims[1]+(xlims[1]*0.03))
            ax.set_title(F'Pathway analysis with KEGG: {title_type}')

        return fig


    def save(self, dir):
        """Save the content of the report as plots.

        :param dir: Path to a directory to save the output to.
        :type dir: string
        """

        # make sure given directory path ends with '/'
        if not dir.endswith('/'):
            dir = dir + '/'

        # collect all produced file in one directory
        try:
            Path(F"{dir}pathway-analysis/").mkdir(parents=True, exist_ok=False)
            print(F'Creating new directory {F"{dir}pathway-analysis/"}')
        except FileExistsError:
            print('Given directory already has required structure.')

        # create and save plots
        # a) for the counts
        if self.total_reac and self.kegg_count:
            count_fig = self.visualise_kegg_counts()
            count_fig.savefig(F'{dir}pathway-analysis/kegg_anno_counts.png', bbox_inches='tight')
        # b) for the actual pathways
        # 1.) global KEGG IDs
        if self.kegg_global:
            # with id
            fig = self.visualise_kegg_pathway(plot_type='global', label='id')
            fig.savefig(F'{dir}pathway-analysis/pathway_global_id.png', bbox_inches='tight')
            # with name
            fig = self.visualise_kegg_pathway(plot_type='global', label='name')
            fig.savefig(F'{dir}pathway-analysis/pathway_global_name.png', bbox_inches='tight')
        # 2.) Overview KEGG IDs
        if self.kegg_over:
            # with id
            fig = self.visualise_kegg_pathway(plot_type='overview', label='id')
            fig.savefig(F'{dir}pathway-analysis/pathway_overview_id.png', bbox_inches='tight')
            # with name
            fig = self.visualise_kegg_pathway(plot_type='overview', label='name')
            fig.savefig(F'{dir}pathway-analysis/pathway_overview_name.png', bbox_inches='tight')
        # 3.) rest
        if self.kegg_paths:
            # grouped by high-level terms
            fig = self.visualise_kegg_pathway(plot_type='high', label='name')
            fig.savefig(F'{dir}pathway-analysis/pathway_high.png', bbox_inches='tight')
            # all with id
            fig = self.visualise_kegg_pathway(plot_type='existing', label='id')
            fig.savefig(F'{dir}pathway-analysis/pathway_existing_id.png', bbox_inches='tight')
            # all with name
            fig = self.visualise_kegg_pathway(plot_type='existing', label='name')
            fig.savefig(F'{dir}pathway-analysis/pathway_existing_name.png', bbox_inches='tight')
