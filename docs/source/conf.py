# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

# Imports
import os
import sys

# Add the source code directory to sys.path
sys.path.insert(0, os.path.abspath('../src'))

# Define project information
project = 'pyissm'
copyright = '2026, ACCESS-NRI'
author = 'ACCESS-NRI'
release = 'v0.0.1.dev190126'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.autosummary',
    'sphinx.ext.napoleon',
    'sphinx_autodoc_typehints',
]

autosummary_generate = True

templates_path = ['_templates']
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']

## TEMPORARY FIX FOR ILL-FORMATTED DOCSTRINGS
# Treat docstring errors as warnings, not fatal errors and do not fail build on docstring formatting
nitpicky = False
suppress_warnings = ["autodoc.import_object", "autodoc", "docutils"]

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'sphinx_rtd_theme'
html_static_path = ['_static']
