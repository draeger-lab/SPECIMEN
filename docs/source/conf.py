# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'SPECIMEN'
copyright = '2024, Carolin Brune and Gwendolyn O. Döbel'
author = 'Carolin Brune and Gwendolyn O. Döbel'
release = '0.0.dev3'

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'src'))

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
              'sphinx.ext.duration',
              'sphinx_rtd_theme',
              'sphinx.ext.autodoc',
              'sphinx.ext.autosectionlabel',
              'sphinx.ext.mathjax',
              'sphinx_copybutton',
              'nbsphinx',
              'IPython.sphinxext.ipython_console_highlighting',
              'sphinxcontrib.bibtex',
              'sphinx.ext.viewcode',
              'sphinxcontrib.icon'
              ]

# For copy buttons in code blocks
copybutton_selector =  "div.copyable pre"
# For citations
bibtex_bibfiles = ['library.bib']

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']
# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
exclude_patterns = ['_build']


# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = 'sphinx_rtd_theme'

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static']

# Change colours in theme for navigation
html_css_files = [
        'https://cdnjs.cloudflare.com/ajax/libs/font-awesome/6.3.0/css/all.min.css',
        'custom_theme.css',
    ]

# Adds logo to documentation page
html_logo = str(Path('images', 'logos', 'specimen-logo-transparent.png'))
html_theme_options = {
    'logo_only': True
}

#Adds logo as favicon to tab
html_favicon = str(Path('images', 'logos', 'specimen-icon-transparent.png'))

# Changes code highlighting
pygments_style = 'blinds-light'

# Make figures numbered
numfig = True

# Explicitly assign the master document
master_doc = 'index'

# -- Autodoc -----------------------------------------------------------------

autodoc_preserve_defaults = True

# we need those to display the code comments otherwise the functions cannot be imported
autodoc_mock_imports = ["psycopg2", 
                        "gffutils",
                        "cplex.exceptions",
                        "cplex",
                        "cobra",
                        "pandas",
                        "libsbml",
                        "numpy",
                        "bioservices",
                        "bioregistry",
                        "bs4",
                        "memote",
                        "tqdm",
                        "psycopg2",
                        "Bio",
                        "sqlalchemy",
                        "ratelimit",
                        "libchebipy",
                        "ols_client",
                        "charges",
                        "click",
                        "cloup",
                        "databases",
                        "yaml",
                        "sortedcontainers",
                        "colorama",
                        "matplotlib",
                        "seaborn",
                        "venn",
                        'refinegems',
                        'xmltodict'
                        ]

