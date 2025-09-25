"""Classes to generate, handle, manipulate and save reports."""

__author__ = "Carolin Brune"

################################################################################
# requirements
################################################################################

import cobra
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import warnings

from pathlib import Path

from refinegems.classes.reports import ModelInfoReport

################################################################################
# variables
################################################################################

################################################################################
# classes
################################################################################


class SpecimenModelInfoReport(ModelInfoReport):
    """A SPECIMEN-specific version of the
    ModelInfoReport for a given model.
    Since this report is a child class of the ModelInfoReport of refineGEMs, 
    it inherits its parameters. Attributes below are the newly addes ones.

    Attributes:
        reac_origin_counts (dict):
            Dictionary with the label of th origin of the reactions as key and the 
            corresponding counts as values.
    """

    def __init__(self, model: cobra.Model) -> None:

        # call the superclass
        super().__init__(model)

        # find out the origin of the reactions
        reac_origin_counts = {
            "via template": 0,
            "refineGEMs based on MetaNetX": 0,
            "refineGEMs based on KEGG": 0,
            "refineGEMs based on BiGG": 0,
            "else": 0,
        }
        for reac in model.reactions:
            # get origin of reaction (based on workflow notation)
            if "creation" in reac.notes.keys():
                reac_origin_counts["via template"] += 1
            elif "created with" in reac.notes.keys():
                reac_origin_counts[reac.notes["created with"]] += 1
            else:
                reac_origin_counts["else"] += 1

        # add new attribute
        self.reac_origin_c = reac_origin_counts

    # extend format table function from parent class
    def format_table(self) -> pd.DataFrame:
        """Extent the functin format_table to include the reaction origin as
        set by SPECIMEN.

        Returns:
            pd.DataFrame:
                The information in table format.
        """
        table = super().format_table()
        table["#reaction origin"] = (
            str(self.reac_origin_c)
            .replace("{", r"")
            .replace("}", r"")
            .replace("'", r"")
        )
        return table

    def visualise(self, color_palette: str = "YlGn") -> tuple[matplotlib.figure.Figure]:
        """Extend the visualisation function to include a graph for the creation type.

        Args:
            - color_palette (str, optional):
                Color palette to be used.
                Defaults to 'YlGn'.

        Returns:
            tuple: 
                Two graphics (1) and (2):
                
                (1) matplotlib.figure.Figure: The original report figure.
                (2) matplotlib.figure.Figure: Report for the creation origin.
        """

        def plot_origin(data, color_palette):

            # create colour gradient
            try:
                cmap = matplotlib.colormaps[color_palette]
            except ValueError:
                warnings.warn('Unknown color palette, setting it to "YlGn"')
                cmap = matplotlib.colormaps["YlGn"]

            # generate the stacked bar plot
            fig, ax = plt.subplots()

            tdata = {}
            for label, count in data.items():
                tdata[label] = np.array([count])

            bottom = np.zeros(1)
            c = 0.2

            for label, count in tdata.items():
                p = ax.barh(["reacs"], count, label=label, left=0.0, color=cmap(c))
                ax.bar_label(p, count, rotation=270)
                bottom += count
                c += 0.2

            plt.legend()

            return fig

        fig1 = super().visualise(color_palette)
        fig2 = plot_origin(self.reac_origin_c, color_palette)

        return (fig1, fig2)

    def save(self, dir: str, color_palette: str = "YlGn") -> None:
        """Save the report and the

        Args:
            - dir (str):
                Path to a directory to save the output files to.
            - color_palette (str, optional):
                Name of a matplotlib colour palette.
                Used as the input for the figures.
                Defaults to 'YlGn'.
        """

        # save the statistics report
        self.format_table().to_csv(Path(dir, f"{self.name}_report.csv"), sep=";")
        # save the visualisation
        figs = self.visualise(color_palette)
        figs[0].savefig(Path(dir, "info_report_vis.png"), bbox_inches="tight", dpi=400)
        figs[1].savefig(Path(dir, "info_origin_vis.png"), bbox_inches="tight", dpi=400)


class DIAMONDReport():
    """ A report for the BLAST runs with DIAMOND.
    Possible for one to three runs.
    
    Attributes:
        - statistics (dict):
            Dictionary with the taxonomy ids as keys and another dictionary
            with counts for 'mapped', 'not best hits' and 'not found' as values. 
        - duplicates (dict):
            Dictionary with the taxonomy ids as keys and a list of duplicate BLAST
            hits as values.
        - below_cutoff (dict):
            Dictionary with the taxonomy ids as keys and a list of BLAST hits with
            an identity below the cutoff as values.
    """
    
    def __init__(
        self,
        statistics: dict,
        duplicates: dict[str, list],
        below_cutoff: dict[str, list]
    ) -> None:
        self.statistics = statistics
        self.duplicates = duplicates
        self.below_cutoff = below_cutoff

    def visualise(self, color_palette: str = "YlGn") -> matplotlib.figure.Figure:
        """Visualise the contents of the report.

        Args:
            - color_palette (str, optional):
                A colour gradient from the matplotlib library.
                If the name does not exist, uses the default.
                Defaults to 'YlGn'.

        Returns:
            If plotting possible: matplotlib.figure.Figure:
                    The plotted figure.
            Else None
        """
        
        try:
            cmap = matplotlib.colormaps[color_palette]
        except ValueError:
            warnings.warn('Unknown color palette, setting it to "YlGn"')
            cmap = matplotlib.colormaps["YlGn"]

        taxids = list(self.statistics.keys())

        fig, ax = plt.subplots()
        positions = np.arange(len(self.statistics[taxids[0]].keys()))

        # Visualise three DIAMOND runs (with three different taxonomy specifications)
        if len(self.statistics.keys()) == 3:
            barWidth = 0.33
            
            bars1 = ax.bar(positions, self.statistics[taxids[0]].values(), label=list(self.statistics.keys())[0], width=0.3, color=cmap(0.25))
            ax.bar_label(bars1, self.statistics[taxids[0]].values())
            bars2 = ax.bar(positions+barWidth, self.statistics[taxids[1]].values(), label=list(self.statistics.keys())[1], width=0.3, color=cmap(0.5)) 
            ax.bar_label(bars2, self.statistics[taxids[1]].values())
            bars3 = ax.bar(positions+2*barWidth, self.statistics[taxids[2]].values(), label=list(self.statistics.keys())[2], width=0.3, color=cmap(0.75))
            ax.bar_label(bars3, self.statistics[taxids[2]].values())
            
            ax.set_ylabel('count', labelpad=10)
            plt.legend()
            plt.xticks(positions+barWidth, ['mapped', 'not best hits', 'not found'])
            plt.title('Report for the DIAMOND runs')

        # Visualise two DIAMOND runs (with two different taxonomy specifications)
        elif len(self.statistics.keys()) == 2:
            barWidth = 0.45
            
            bars1 = ax.bar(positions, self.statistics[taxids[0]].values(), label=list(self.statistics.keys())[0], width=0.4, color=cmap(0.5)) 
            ax.bar_label(bars1, self.statistics[taxids[0]].values())
            bars2 = ax.bar(positions+barWidth, self.statistics[taxids[1]].values(), label=list(self.statistics.keys())[1], width=0.4, color=cmap(0.75))
            ax.bar_label(bars2, self.statistics[taxids[1]].values())
            
            ax.set_ylabel('count', labelpad=10)
            plt.legend()
            plt.xticks(positions+barWidth/2, ['mapped', 'not best hits', 'not found'])
            plt.title('Report for the DIAMOND runs')
        
        # Visualise one DIAMOND run
        elif len(self.statistics.keys()) == 1:
            bars = ax.bar(self.statistics[taxids[0]].keys(), self.statistics[taxids[0]].values(), label=list(self.statistics.keys())[0], width=0.7, color=cmap(0.75))
            ax.bar_label(bars, self.statistics[taxids[0]].values())
            
            ax.set_ylabel('count', labelpad=10)
            plt.legend()
            plt.xticks(positions, ['mapped', 'not best hits', 'not found'])
            plt.title('Report for the DIAMOND run')
        
        return fig

    def save(self, dir: str, color_palette: str = 'YlGn'):
        """Save the report in a new folder called 'DIAMOND_report'.
        A visualization of the counts as well as tsv files for duplicates
        and hits below the cutoff for each of the runs are saved.

        Args:
            - dir (str):
                Path to a directory to save the report to.
            - color_palette (str, optional):
                A colour gradient from the matplotlib library.
                If the name does not exist, uses the default.
                Defaults to 'YlGn'.
        """
        
        dir = str(Path(dir, 'DIAMOND_report'))
        Path(dir).mkdir(parents=True, exist_ok=False)

        fig = self.visualise(color_palette)
        fig.savefig(str(Path(dir, 'DIAMOND_visual.png')))
        
        for key in self.duplicates:
            pd.DataFrame(self.duplicates[key]).to_csv(str(Path(dir,f'duplicates_{key}.tsv')), sep="\t")
        for key in self.below_cutoff:
            pd.DataFrame(self.below_cutoff[key]).to_csv(str(Path(dir,f'below_cutoff_{key}.tsv')), sep="\t")