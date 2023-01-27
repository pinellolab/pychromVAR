import inspect
import logging
import os
import sys
from datetime import datetime
from pathlib import Path

from sphinx.application import Sphinx
from sphinx.ext import autosummary

sys.path.insert(0, os.path.abspath('.'))

# -- Project information

project = 'pychromVAR'
author = 'Zhijian Li'

# -- General configuration

extensions = [
    'sphinx.ext.duration',
    'sphinx.ext.doctest',
    'sphinx.ext.autodoc',
    'sphinx.ext.autosummary',
    'sphinx.ext.intersphinx',
    'sphinx.ext.autosectionlabel',
    'numpydoc',
    'nbsphinx',
    'IPython.sphinxext.ipython_console_highlighting'
]

intersphinx_mapping = {
    'python': ('https://docs.python.org/3/', None),
    'sphinx': ('https://www.sphinx-doc.org/en/master/', None),
}
intersphinx_disabled_domains = ['std']

templates_path = ['_templates']

# -- Options for HTML output
master_doc = 'index'

html_theme = 'sphinx_rtd_theme'
html_static_path = ["_static"]
html_theme_options = dict(
    logo_only=True,
    display_version=True,
)
html_context = dict(
    display_github=True,  # Integrate GitHub
    github_user='saezlab',  # Username
    github_repo='decoupler-py',  # Repo name
    github_version='master',  # Version
    conf_py_path='/docs/source/',  # Path in the checkout to the docs root
)
html_show_sphinx = False
html_logo = 'logo.png'
html_favicon = 'logo.png'
html_css_files = [
    'css/custom.css',
]

# -- Options for EPUB output
epub_show_urls = 'footnote'