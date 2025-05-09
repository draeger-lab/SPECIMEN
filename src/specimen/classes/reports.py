"""Classes to generate, handle, manipulate and save reports."""

__author__ = "Carolin Brune"

################################################################################
# requirements
################################################################################

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

    Attributes:
        model:
            The GEM loaded with COBRApy.
    """

    def __init__(self, model) -> None:

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
