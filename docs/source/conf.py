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
    "breathe",
    "myst_parser",
    "sphinx.ext.intersphinx",
    "sphinx_sitemap",
    "sphinx_design",
    "sphinxcontrib.mermaid",
]

breathe_projects = {"seams": "xml"}
breathe_default_project = "seams"

templates_path = ["_templates"]
exclude_patterns = []

intersphinx_mapping = {
    "python": ("https://docs.python.org/3", None),
    "numpy": ("https://numpy.org/doc/stable", None),
    "ase": ("https://wiki.fysik.dtu.dk/ase/", None),
    "pydseams": ("https://d-seams.github.io/PydSEAMSlib/", None),
    "luadseams": ("https://d-seams.github.io/yodaStruct/", None),
}

# -- Options for HTML output -------------------------------------------------
html_theme = "shibuya"
html_static_path = ["_static"]
html_css_files = ["custom.css"]
html_js_files = []
html_favicon = "_static/dSeamsLogo.ico"

# Mermaid: default CDN; diagrams via ``.. mermaid::`` from Org RST export.
mermaid_version = "11.4.0"
mermaid_init_js = "mermaid.initialize({startOnLoad:true, theme:'neutral'});"

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
    "light_logo": "_static/dSeamsLogo.png",
    "dark_logo": "_static/dSeamsLogo.png",
    "globaltoc_expand_depth": 1,
    "toctree_collapse": True,
    "toctree_maxdepth": 3,
    "toctree_titles_only": True,
    "nav_links": [
        {
            "title": "Ecosystem",
            "children": [
                {
                    "title": "d-SEAMS engine",
                    "url": "https://docs.dseams.info",
                    "summary": "libyodaLib and the seams CLI",
                    "external": True,
                },
                {
                    "title": "pydseams",
                    "url": "https://d-seams.github.io/PydSEAMSlib/",
                    "summary": "Python Frame API on yoda",
                    "external": True,
                },
                {
                    "title": "dseams (Lua)",
                    "url": "https://d-seams.github.io/yodaStruct/",
                    "summary": "require(\"dseams\") and Fennel",
                    "external": True,
                },
            ],
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

html_baseurl = "https://docs.dseams.info/"
