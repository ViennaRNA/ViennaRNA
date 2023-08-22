# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
import os
import sys
sys.path.insert(0, os.path.abspath('.'))
sys.path.insert(0, os.path.abspath('../'))
sys.path.insert(0, os.path.abspath('../../'))


# -- Project information -----------------------------------------------------

project = 'RNA'
copyright = '2022, Ronny Lorenz'
author = 'Ronny Lorenz'

# The full version, including alpha/beta/rc tags
release = '2.5.0'


# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
#'numpydoc'
    'sphinx.ext.napoleon',
    'sphinx.ext.autodoc',
    'sphinx.ext.viewcode',
    'sphinx.ext.autosummary',
    'sphinx.ext.graphviz',
    'sphinx_rtd_theme',
    'myst_parser',
    'sphinx.ext.mathjax'
]

# Sort members by type
#autodoc_member_order = 'groupwise'

napoleon_google_docstring = False
napoleon_numpy_docstring = True
napoleon_include_init_with_doc = False
napoleon_include_private_with_doc = False
napoleon_include_special_with_doc = True
napoleon_use_admonition_for_examples = False
napoleon_use_admonition_for_notes = False
napoleon_use_admonition_for_references = False
napoleon_use_ivar = False
napoleon_use_param = True
napoleon_use_rtype = True
napoleon_preprocess_types = True
napoleon_type_aliases = None
napoleon_attr_annotations = True
mathjax_path = "js/mathjax/tex-chtml.js"
mathjax3_config = {
    'tex': {
        'inlineMath': [['$', '$'], ['\\(', '\\)']],
        'displayMath': [["\\[", "\\]"]]
    }
}
napoleon_type_aliases = {
    "PRIVATE int": "int",
    "PRIVATE FLT_OR_DBL" : "double",
    "unsigned int *": "list-like(unsigned int)",
    "unsigned int **": "list-like(list-like(unsigned int))",
    "short *": "list-like(int)",
    "char *": "string",
    "const char *": "string",
    "float *": "list-like(double)",
    "double *": "list-like(double)",
    "double **": "list-like(list-like(double))",
    "vrna_fold_compound_t *" : "fold_compound",
    "vrna_param_t *" : "param",
    "vrna_exp_param_t *" : "exp_param",
    "vrna_md_t *" : "md",
    "std::string": "string",
    "FLT_OR_DBL" : "double",
    "FLT_OR_DBL *" : "list-like(double)",
    "std::vector<FLT_OR_DBL>" : "list-like(double)"
}

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# The encoding of source files.
source_encoding = "utf-8"

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = []


# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = 'sphinx_rtd_theme'

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static']
