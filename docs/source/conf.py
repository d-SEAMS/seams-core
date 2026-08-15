#!/usr/bin/env python3

import os
import sys

# -- Path setup --------------------------------------------------------------
sys.path.insert(0, os.path.abspath("../../subprojects/doxyrest/sphinx"))

# -- Project information -----------------------------------------------------
project = "d-SEAMS"
copyright = "2019-2026, d-SEAMS core team"
author = "Rohit Goswami, Amrita Goswami, Ruhila S"

# -- General configuration ---------------------------------------------------
extensions = [
    "doxyrest",
    "cpplexer",
    "myst_parser",
    "sphinx.ext.intersphinx",
    "sphinx_sitemap",
]

templates_path = ["_templates"]
exclude_patterns = []

intersphinx_mapping = {
    "python": ("https://docs.python.org/3", None),
    "numpy": ("https://numpy.org/doc/stable", None),
    "ase": ("https://wiki.fysik.dtu.dk/ase/", None),
    "pydseams": ("https://d-seams.github.io/PydSEAMSlib/", None),
}

# -- Options for HTML output -------------------------------------------------
html_theme = "shibuya"
html_static_path = ["_static"]

html_context = {
    "source_type": "github",
    "source_user": "d-SEAMS",
    "source_repo": "seams-core",
    "source_version": "main",
    "source_docs_path": "/docs/source/",
}

html_theme_options = {
    "github_url": "https://github.com/d-SEAMS/seams-core",
    "accent_color": "teal",
    "dark_code": True,
    "globaltoc_expand_depth": 1,
    "nav_links": [
        {
            "title": "pydseams",
            "url": "https://github.com/d-SEAMS/PydSEAMSlib",
            "external": True,
        },
        {
            "title": "dseams (Lua)",
            "url": "https://github.com/d-SEAMS/yodaStruct",
            "external": True,
        },
    ],
}

html_sidebars = {
    "**": [
        "sidebars/localtoc.html",
        "sidebars/repo-stats.html",
        "sidebars/edit-this-page.html",
    ],
}

html_baseurl = "https://docs.dseams.info"
