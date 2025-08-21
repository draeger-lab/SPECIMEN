[![GitHub License](https://img.shields.io/github/license/draeger-lab/specimen)](https://opensource.org/license/GPL-3.0)
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)
![Python Version from PEP 621 TOML](https://img.shields.io/python/required-version-toml?tomlFilePath=https%3A%2F%2Fraw.githubusercontent.com%2Fdraeger-lab%2Fspecimen%2Fmain%2Fpyproject.toml)
[![Documentation Status](https://readthedocs.org/projects/specimen/badge/?version=latest)](https://specimen.readthedocs.io/en/latest/?badge=latest)
![GitHub release (with filter)](https://img.shields.io/github/v/release/draeger-lab/specimen?logo=github&label=SPECIMEN&color=B4A069&style=flat-square&include_prereleases)
![GitHub last commit (branch)](https://img.shields.io/github/last-commit/draeger-lab/specimen/main)
![Repo Size](https://img.shields.io/github/repo-size/draeger-lab/specimen)
![GitHub all releases](https://img.shields.io/github/downloads/draeger-lab/specimen/total?logo=github&label=GitHub%20downloads)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.12723500.svg)](https://doi.org/10.5281/zenodo.12723500)

![Logo of SPECIMEN](docs/source/images/LogoSPECIMEN.png)

# SPECIMEN

SPECIMEN is a collection of different workflows designed for the automated and standardised curation of genome-scale models. It is mainly based on the [refineGEMs toolbox](https://github.com/draeger-lab/refinegems/tree/main), but also includes additional tools like [CarveMe](https://carveme.readthedocs.io/en/latest/).

> Note: <br> 
  This tool is currently still under active developement, any feedback or ideas are welcome. Also, since its developemental, please note that bugs and error will probably occur. Please feel free to report them to the developers.

Currently avaible workflow:

- CMPB - CarveMe-ModelPolisher based:\\Starting from a CarveMe draft model, refine and extend it towards a high-quality stain-specific model

- HQTB - High-quality template based:\\This pipeline follows the modelling approach of using a high-quality template model as a basis for the reconstruction of a new model from a new genome (e.g. a different strain). 

> Note: Due to some major refactoring changes in refineGEMs this workflow might not run as expected.

- PGAB: under construction

<!-- TOC -->
- [Installation](#installation)
    - [Via pip](#pypi-via-pip)
    - [Via Docker](#docker-via-docker)
- [Quickstart](#quickstart)
    - [After install via pip](#pypi-after-install-via-pip)
    - [After install via Docker](#docker-after-install-via-docker)
- [Documentation](#documentation)
- [Repositories using SPECIMEN](#repositories-using-specimen)

<!-- /TOC -->

## Installation  
The workflow collection ``SPECIMEN``can be installed via pip or via Docker.

### ![pypi](https://skillicons.dev/icons?i=py) Via pip 

Download this repository and run the command `pip install -e .` inside the top-level directory.     

When running certain steps, further tools need to be installed:

- [DIAMOND, version 2.0.4 or higher](https://github.com/bbuchfink/diamond)
- [EntrezDirect](https://www.ncbi.nlm.nih.gov/books/NBK179288/), if no NCBI mapping has been created beforehand

### ![docker](https://skillicons.dev/icons?i=docker) Via Docker
``SPECIMEN`` can also be used via Docker. To build the docker image, firtsly clone the repository:

```bash
   git clone "https://github.com/draeger-lab/specimen.git"
```

Then change into the directory and build the image:

```bash
   cd specimen \
   docker build -t specimen .
```

## Quickstart

### ![pypi](https://skillicons.dev/icons?i=py) After install via pip

After the installation, main functionalities can be accessed either via the command line. Try running `specimen --help` for more information.

For greater control or for further integration into other scripts, the modules of SPECIMEN can be loaded as a Python package using `import specimen` in a Python script.

### ![docker](https://skillicons.dev/icons?i=docker) After install via Docker
The default command executed by the image is ``specimen -h`` and provides the help information for the CLI of 
``SPECIMEN``.

```bash
   docker run specimen
```

To use the image interactively and open a bash shell, run the following command:

```bash
   docker run -it --entrypoint bash specimen
```

To use the image for specific commands, you can simply use every of the CLI commands as entrypoint. 
For example, to run the CMPB pipeline, use:

```bash
   docker run --name <container_name> -v <user_folder>:/sp_cont specimen cmpb run ./path/to/CMPB_config.yaml
```

## Documentation

> [!WARNING]
> 🚧 The documentation is currently under heavy-rework!

For more information about the available pipelines, the code or for troubleshooting, please refer to the documentation of the tool [here](https://specimen.readthedocs.io/en/latest/).

## Repositories using SPECIMEN
- draeger-lab/Cacnes - `private`
- draeger-lab/Cgranulosum - `private`
- draeger-lab/Koxytoca - `private`
- draeger-lab/Mfortuitum - `private`
- draeger-lab/Scohnii - `private`
