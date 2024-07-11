[![GitHub License](https://img.shields.io/github/license/draeger-lab/specimen)](https://opensource.org/license/GPL-3.0)
![Python Version from PEP 621 TOML](https://img.shields.io/python/required-version-toml?tomlFilePath=https%3A%2F%2Fraw.githubusercontent.com%2Fdraeger-lab%2Fspecimen%2Fmain%2Fpyproject.toml)
[![Documentation Status](https://readthedocs.org/projects/specimen/badge/?version=latest)](https://specimen.readthedocs.io/en/latest/?badge=latest)

![GitHub release (with filter)](https://img.shields.io/github/v/release/draeger-lab/specimen?logo=github&label=SPECIMEN&color=B4A069&style=flat-square&include_prereleases)

![GitHub last commit (branch)](https://img.shields.io/github/last-commit/draeger-lab/specimen/main)
![Repo Size](https://img.shields.io/github/repo-size/draeger-lab/specimen)
![GitHub all releases](https://img.shields.io/github/downloads/draeger-lab/specimen/total?logo=github&label=GitHub%20downloads)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.12723501.svg)](https://doi.org/10.5281/zenodo.12723501)

![Logo of SPECIMEN](docs/source/images/LogoSPECIMEN.png)

# SPECIMEN

SPECIMEN is a collection of different workflows designed for the automated and standardised curation of genome-scale models. It is mainly based on the [refineGEMs toolbox](https://github.com/draeger-lab/refinegems/tree/dev-2), but also includes additional tools like [CarveMe](https://carveme.readthedocs.io/en/latest/).

> Note: This tool is currently still under active developement, any feedback or ides are welcome.

Currently avaible workflow:

- CMPB - CarveMe-ModelPolisher based:\\Starting from a CarveMe draft model, refine and extend it towards a high-quality stain-specific model

    > Note: For a future update, optional direct integration of CarveMe into the pipeline is planned

- HQTB - High-quality template based:\\This pipeline follows the modelling approach of using a high-quality template model as a basis for the reconstruction of a new model from a new genome (e.g. a different strain). 

- PGAB: under construction

## Installation  

Download this repository and run the command `pip install -e .` inside the top-level directory.     

When running the HQTB Pipeline, further tools need to be installed:

- [DIAMOND, version 2.0.4 or higher](https://github.com/bbuchfink/diamond)
- [EntrezDirect](https://www.ncbi.nlm.nih.gov/books/NBK179288/), if no NCBI mapping has been created beforehand

## Quickstart

After the installation, main functionalities can be accessed either via the command line. Try running `specimen --help` for more information.

For greater control or for further integration into other scripts, the modules of SPECIMEN can be loaded as a Python package using `import specimen` in a Python script.

## Documentation

For more information about the available pipelines, the code or for troubleshooting, please refer to the documentation of the tool [here](https://specimen.readthedocs.io/en/latest/).
