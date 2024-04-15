# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'SPECIMEN'
copyright = '2024, Carolin Brune'
author = 'Carolin Brune'
release = '1.0.0'

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
import os
import sys
sys.path.insert(0, os.path.abspath('../../src/'))

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
              'sphinx.ext.duration',
              'sphinx_rtd_theme',
              'sphinx.ext.autodoc',
              'sphinx.ext.autosectionlabel',
              'sphinx.ext.mathjax',
              'sphinx_copybutton',
              'nbsphinx',
              'IPython.sphinxext.ipython_console_highlighting',
              'sphinxcontrib.bibtex'
              ]

templates_path = ['_templates']
exclude_patterns = []

# For copy buttons in code blocks
copybutton_selector =  "div.copyable pre"
# For citations
bibtex_bibfiles = ['library.bib']


# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'sphinx_rtd_theme'
html_static_path = ['_static']

# Change colours in theme for navigation
html_css_files = ['custom_theme.css']

# Adds logo to documentation page
html_logo = 'images/LogoSPECIMEN.png'
html_theme_options = {
    'logo_only': True,
    'display_version': False
}

#Adds logo as favicon to tab
html_favicon = 'images/LogoSPECIMEN.png'

# Changes code highlighting
pygments_style = 'blinds-light'

# Make figures numbered
numfig = True

# Explicitly assign the master document
master_doc = 'index'