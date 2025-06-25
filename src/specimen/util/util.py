#!/usr/bin/env python
"""Utility functions."""

__author__ = "Carolin Brune"

################################################################################
# requirements
################################################################################

import os
import os.path
import subprocess
import sys

from Bio import SeqIO
from pathlib import Path

from refinegems.utility.io import parse_gbff_for_cds


################################################################################
# variables
################################################################################

################################################################################
# functions
################################################################################

# create a DIAMOND reference database
# -----------------------------------


def create_DIAMOND_db_from_folder(
    dir: str, out: str, name: str = "database", extension: str = "faa", threads: int = 2
):
    """Build a DIAMOND database from a folder containing FASTA files.

    Args:
        - dir (str):
            Path to the directory to search for FASTA files for
            the database (recursive file search).
        - out (str):
            Path of the directory of the output.
        - name (str, optional):
            Name of the created database.
            Defaults to 'database'.
        - extension (str, optional):
            File extension of the FASTA files
            (to determine which files to search for).
            Defaults to 'faa'.
        - threads (int, optional):
            Number of threads to use for DIAMOND.
            Defaults to 2.
    """

    # get fasta file names
    fasta_files = Path(dir).rglob(f"*.{extension}")

    # -----------------------
    # combine to 1 FASTA file
    # -----------------------

    # check if folder already has a combined FASTA
    outname_fasta = Path(out, "combinded.faa")
    save = True
    if os.path.isfile(outname_fasta):
        print("A combined.faa files already exists in the given folder.")
        answer = None
        while answer not in ("overwrite", "use", "end"):
            answer = input(
                "Do you want to [overwrite] or [use] the file or [end] the programm?: "
            )
            if answer == "use":
                save = False
            elif answer == "end":
                sys.exit()
            elif answer == "overwrite":
                save = True
            else:
                print("Enter [overwrite] or [use] or [end] (without brackets)")

    # save all in new (combined) FASTA
    if save:
        with open(outname_fasta, "w") as outfile:
            for f in fasta_files:
                SeqIO.write(SeqIO.parse(f, "fasta"), outfile, "fasta")

    # -------------------------
    # generate DIAMOND database
    # -------------------------

    outname_dnmd = Path(out, name)
    bl = "\\ "
    print(
        f'running the following command:\ndiamond makedb --in {str(outname_fasta).replace(" ",bl)} -d {str(outname_dnmd).replace(" ",bl)} -p {threads}'
    )
    subprocess.run(
        [
            f'diamond makedb --in {str(outname_fasta).replace(" ",bl)} -d {str(outname_dnmd).replace(" ",bl)} -p {threads}'
        ],
        shell=True,
    )
    print("finished")


# create a NCBI mapping file for the database
# -------------------------------------------


def create_NCBIinfo_mapping(dir: str, out: str, extension: Literal["gbff"] = "gbff"):
    """Create a NCBI information mapping file from a folder containing e.g. gbff files.

    Args:
        - dir (str):
            Path to the directory for the recursive file search for the mapping.
        - out (str):
            Path of the directory for the output.
        - extension (Literal['gbff'], optional):
            Name of the file extension to be searched.
            Default is gbff, and currently it is advised to leave it at that.
            Defaults to 'gbff'.
    """

    # get gbff file names
    list_files = Path(dir).rglob(f"*.{extension}")
    file_no = len(list(Path(dir).rglob(f"*.{extension}")))

    # -----------------------------------------------
    # extract information and create the mapping file
    # -----------------------------------------------

    with open(out, "w") as out_file:

        file_counter = 1

        for file_name in list_files:

            print(f"Parsing {file_counter}/{file_no}: {file_name}")

            # retrieve information
            match extension:
                case "gbff":
                    info = parse_gbff_for_cds(file_name)
                case "gff":
                    raise NotImplementedError(
                        "This option is under discussion and has yet to be implemented."
                    )

                case _:
                    raise ValueError(f"Unknown file extension {extension}")

            # save
            info.to_csv(out_file, header=True, index=False)
            # go to next
            file_counter += 1
