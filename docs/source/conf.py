
# Importando diretório
import os
import sys
sys.path.insert(0, os.path.abspath('../../'))

# sphinx-build -M html docs/source/ docs/build/

# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'dopings'
copyright = '2024, Gustavo Campos'
author = 'Gustavo Campos'
release = '0.1'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = ['sphinx.ext.autodoc', 'sphinx.ext.doctest', 
              'sphinx.ext.todo', 'sphinx.ext.ifconfig', 
              'sphinx.ext.githubpages', 'sphinx.ext.napoleon' ]

templates_path = ['_templates']
exclude_patterns = []

language = 'en'

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = "sphinx_rtd_theme"
html_static_path = ['_static']
