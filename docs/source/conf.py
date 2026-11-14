#!/usr/bin/env python3

import os
import sys

# -- Path setup --------------------------------------------------------------
_source = os.path.dirname(os.path.abspath(__file__))
_docs = os.path.dirname(_source)


def _write_sitting_tables():
    """Fill the performance tables from the CSV before Sphinx reads pages."""
    import importlib.util

    script = os.path.join(_docs, "render_sitting.py")
    spec = importlib.util.spec_from_file_location("render_sitting", script)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    if module.main([]) != 0:
        raise SystemExit("rc3 sitting tables do not match the CSV")


_write_sitting_tables()


def _write_h2h():
    """Fill the performance head-to-head from the benchmark CSVs."""
    import importlib.util

    script = os.path.join(_docs, "render_h2h.py")
    spec = importlib.util.spec_from_file_location("render_h2h", script)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    if module.main(["--check"]) != 0:
        raise SystemExit("performance head-to-head does not match the CSVs")


_write_h2h()


def _write_comparisons():
    """Fill the Comparisons tables from the vendored benchmark CSVs."""
    import importlib.util

    script = os.path.join(_docs, "render_comparisons.py")
    spec = importlib.util.spec_from_file_location("render_comparisons", script)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    if module.main(["--check"]) != 0:
        raise SystemExit("Comparisons tables do not match the vendored CSVs")


_write_comparisons()


def _write_performance_range():
    """Fill the 4096 to 110592 neighbour and ring tables from the CSVs."""
    import importlib.util

    script = os.path.join(_docs, "render_performance_range.py")
    spec = importlib.util.spec_from_file_location("render_performance_range", script)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    if module.main(["--check"]) != 0:
        raise SystemExit("performance range does not match the CSVs")


_write_performance_range()


def _write_science_gallery():
    """Fill the four-figure page from the shared census CSV."""
    import importlib.util

    script = os.path.join(_docs, "write_science_gallery.py")
    spec = importlib.util.spec_from_file_location("write_science_gallery", script)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    if module.main([]) != 0:
        raise SystemExit("science gallery was not written")


_write_science_gallery()


def _write_gallery():
    """Draw the science gallery from the slow tutorial runs."""
    import importlib.util

    script = os.path.join(_docs, "render_gallery.py")
    spec = importlib.util.spec_from_file_location("render_gallery", script)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    try:
        status = module.main([])
    except SystemExit as exc:
        status = exc.code if isinstance(exc.code, int) else 1
    if status not in (0, None):
        kept = os.path.isdir(os.path.join(_docs, "source", "tutorials", "gallery"))
        if not kept:
            raise SystemExit("science gallery failed")
        print("sphinx: science gallery left as committed", file=sys.stderr)


_write_gallery()


def _write_molecules():
    """Draw chemparseplot stills before Sphinx reads the figure list."""
    import importlib.util

    script = os.path.join(_docs, "scripts", "render_molecules.py")
    spec = importlib.util.spec_from_file_location("render_molecules", script)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    argv = ["--out", os.path.join(_source, "tutorials", "molecules"), "--pages"]
    libdir = os.environ.get("SEAMS_LIBDIR", "")
    if libdir:
        argv.extend(["--libdir", libdir])
    try:
        status = module.main(argv)
    except SystemExit as exc:
        status = exc.code if isinstance(exc.code, int) else 1
    if status not in (0, None):
        if not os.path.isdir(os.path.join(_source, "tutorials", "molecules")):
            raise SystemExit("molecule stills failed")
        print("sphinx: molecule stills left as committed", file=sys.stderr)


_write_molecules()


def _write_video_players():
    """Fill the Odysee player list from docs/videos/<slug>/embed.txt."""
    import importlib.util
    import pathlib
    import shutil

    script = os.path.join(_docs, "video_page.py")
    spec = importlib.util.spec_from_file_location("video_page", script)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    videos = pathlib.Path(_docs) / "videos"
    dest = pathlib.Path(_source) / "videos" / "players.rst"
    module.write_players(videos, dest)
    static = os.path.join(_source, "_static", "captions")
    os.makedirs(static, exist_ok=True)
    for name in os.listdir(videos):
        track = os.path.join(videos, name, "captions.vtt")
        if os.path.isfile(track):
            shutil.copy(track, os.path.join(static, name + ".vtt"))


_write_video_players()


def _write_pool_v5():
    """Fill the pool-audit v5 tables before Sphinx reads the benchmarks page."""
    import importlib.util

    script = os.path.join(_docs, "render_pool_v5.py")
    spec = importlib.util.spec_from_file_location("render_pool_v5", script)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    if module.main([]) != 0:
        raise SystemExit("pool-audit v5 tables do not match the log")


_write_pool_v5()


def _write_release_notes():
    """Fill the 2.14.0 notes from changelog.d before Sphinx reads pages."""
    import importlib.util

    script = os.path.join(_docs, "check_changelog.py")
    spec = importlib.util.spec_from_file_location("check_changelog", script)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    if module.main([]) != 0:
        raise SystemExit("2.14.0 notes do not match changelog.d")


_write_release_notes()


def _write_padding():
    """Fill the cutoff-and-padding how-to from the validation CSV."""
    import importlib.util

    script = os.path.join(_docs, "render_padding.py")
    spec = importlib.util.spec_from_file_location("render_padding", script)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    if module.main([]) != 0:
        raise SystemExit("cutoff and padding tables do not match the CSV")


_write_padding()


def _stage_markdown_book():
    """Copy the language pages into the Sphinx tree so myst renders them."""
    import pathlib

    groups = (
        ("api", "stable", "Stable API"),
        ("engines", "engines", "Engine builds"),
        ("tutorials", "fronts", "Language tutorials"),
    )
    for src_name, dest_name, title in groups:
        src = pathlib.Path(_docs) / src_name
        dest = pathlib.Path(_source) / dest_name
        dest.mkdir(parents=True, exist_ok=True)
        names = []
        for page in sorted(src.glob("*.md")):
            (dest / page.name).write_text(page.read_text(encoding="utf-8"), encoding="utf-8")
            names.append(page.stem)
        underline = "=" * len(title)
        # The sidebar toctree lives on the landing page. A second toctree
        # here would put every page in two trees.
        body = [title, underline, ""]
        for name in names:
            body.append("- :doc:`{}`".format(name))
        body.append("")
        (dest / "index.rst").write_text("\n".join(body), encoding="utf-8")


_stage_markdown_book()


def _link_capi():
    """Point the generated book at the committed C API section."""
    index = os.path.join(_source, "index.rst")
    if not os.path.isfile(index):
        return
    text = open(index, encoding="utf-8").read()
    if "\n   capi/index\n" in text:
        return
    extra = (
        "\n.. toctree::\n"
        "   :hidden:\n"
        "   :maxdepth: 2\n"
        "   :caption: C API\n"
        "\n"
        "   capi/index\n"
    )
    with open(index, "a", encoding="utf-8") as handle:
        handle.write(extra)


_link_capi()


def _link_lang():
    """Point the generated book at the per-language API pages."""
    index = os.path.join(_source, "index.rst")
    if not os.path.isfile(index):
        return
    text = open(index, encoding="utf-8").read()
    if "lang/index" in text:
        return
    extra = (
        "\n.. toctree::\n"
        "   :hidden:\n"
        "   :maxdepth: 2\n"
        "   :caption: Language APIs\n"
        "\n"
        "   lang/index\n"
    )
    with open(index, "a", encoding="utf-8") as handle:
        handle.write(extra)


_link_lang()


def _link_science():
    """Point the generated book at the science tutorials."""
    index = os.path.join(_source, "index.rst")
    if not os.path.isfile(index):
        return
    text = open(index, encoding="utf-8").read()
    if "science/index" in text:
        return
    extra = (
        "\n.. toctree::\n"
        "   :hidden:\n"
        "   :maxdepth: 2\n"
        "   :caption: Science tutorials\n"
        "\n"
        "   science/index\n"
    )
    with open(index, "a", encoding="utf-8") as handle:
        handle.write(extra)


def _science_figure_stubs():
    """Keep the book buildable before the tutorial script has drawn."""
    dest = os.path.join(_source, "science", "figures")
    os.makedirs(dest, exist_ok=True)
    stub = (
        '<svg xmlns="http://www.w3.org/2000/svg" width="640" height="80">'
        '<rect width="100%" height="100%" fill="white"/>'
        '<text x="16" y="44" font-size="14">'
        "Figure is drawn by docs/examples/science_tutorials.py"
        "</text></svg>\n"
    )
    for name in ("brine", "cages"):
        path = os.path.join(dest, name + ".svg")
        if not os.path.isfile(path):
            with open(path, "w", encoding="utf-8") as handle:
                handle.write(stub)
    for name in ("interface", "brine", "sii", "sh"):
        _write_png_stub(os.path.join(dest, name + ".png"))


def _write_png_stub(png):
    if os.path.isfile(png):
        return
    import struct
    import zlib

    raw = b"\x00\xff\xff\xff"

    def chunk(tag, data):
        return (
            struct.pack(">I", len(data))
            + tag
            + data
            + struct.pack(">I", zlib.crc32(tag + data) & 0xFFFFFFFF)
        )

    blob = b"\x89PNG\r\n\x1a\n" + chunk(b"IHDR", struct.pack(">IIBBBBB", 1, 1, 8, 2, 0, 0, 0))
    blob += chunk(b"IDAT", zlib.compress(raw))
    blob += chunk(b"IEND", b"")
    with open(png, "wb") as handle:
        handle.write(blob)


_link_science()
_science_figure_stubs()


def _link_walkthrough():
    """Point the generated book at the engine walkthroughs."""
    index = os.path.join(_source, "index.rst")
    if not os.path.isfile(index):
        return
    text = open(index, encoding="utf-8").read()
    if "walkthrough/index" in text:
        return
    extra = (
        "\n.. toctree::\n"
        "   :hidden:\n"
        "   :maxdepth: 2\n"
        "   :caption: Engine walkthroughs\n"
        "\n"
        "   walkthrough/index\n"
    )
    with open(index, "a", encoding="utf-8") as handle:
        handle.write(extra)


def _walkthrough_stubs():
    dest = os.path.join(_source, "walkthrough", "figures")
    os.makedirs(dest, exist_ok=True)
    stub = (
        '<svg xmlns="http://www.w3.org/2000/svg" width="420" height="80">'
        '<rect width="100%" height="100%" fill="white"/>'
        '<text x="16" y="44" font-size="14">engine census</text></svg>\n'
    )
    for name in ("lammps", "gromacs", "openmm", "plumed"):
        path = os.path.join(dest, name + ".svg")
        if not os.path.isfile(path):
            with open(path, "w", encoding="utf-8") as handle:
                handle.write(stub)


_link_walkthrough()
_walkthrough_stubs()
os.environ["DSEAMS_DOC_BUILD"] = "1"
sys.path.insert(0, os.path.abspath(os.path.join(_docs, "..", "bindings", "python")))
_doxyrest = os.path.abspath(os.path.join(_source, "../../subprojects/doxyrest/sphinx"))
_xml_index = os.path.join(_source, "xml", "index.xml")
_api_index = os.path.join(_source, "api", "index.rst")
# Doxygen output is not part of the tree. A missing C++ index must not
# fail the book: write a short page so the Reference toctree resolves.
if not os.path.isfile(_api_index):
    os.makedirs(os.path.dirname(_api_index), exist_ok=True)
    with open(_api_index, "w", encoding="utf-8") as _api_page:
        _api_page.write(
            "API Reference\n"
            "=============\n"
            "\n"
            "These pages are filled by Doxygen and Doxyrest.\n"
            "Until that step has run, the header map is the API page.\n"
        )

# -- Project information -----------------------------------------------------
project = "d-SEAMS"
copyright = "2019-2026, d-SEAMS core team"
author = "Rohit Goswami, Amrita Goswami, Ruhila S"

# -- General configuration ---------------------------------------------------
extensions = [
    "myst_parser",
    "sphinx.ext.autodoc",
    "sphinx.ext.doctest",
    "sphinx.ext.intersphinx",
    "numpydoc",
    "sphinx_sitemap",
    "sphinx_design",
    "sphinxcontrib.mermaid",
]
numpydoc_show_class_members = False
numpydoc_class_members_toctree = False
numpydoc_validation_checks = set()
doctest_global_setup = "import sys\nsys.path.insert(0, {!r})\n".format(
    os.path.join(_docs, "examples")
)
# doxyrest and the XML tree are produced by the docs task. A build from
# the committed sources does not require them to already exist.
if os.path.isfile(_xml_index) and os.path.isdir(_doxyrest):
    sys.path.insert(0, _doxyrest)
    extensions = ["doxyrest", "cpplexer", "breathe", *extensions]
    breathe_projects = {"seams": os.path.join(_source, "xml")}
    breathe_default_project = "seams"

templates_path = ["_templates"]
exclude_patterns = ["science/figures/**"]
# A warning is a failed book. Sphinx stops instead of publishing it.
warningiserror = True

intersphinx_mapping = {
    "python": ("https://docs.python.org/3", None),
    "numpy": ("https://numpy.org/doc/stable", None),
    "pydseams": ("https://d-seams.github.io/PydSEAMSlib/", None),
    "luadseams": ("https://d-seams.github.io/yodaStruct/", None),
}

# -- Options for HTML output -------------------------------------------------
html_theme = "shibuya"
# doi.org is the citation link. These publishers answer the redirect with
# 403 to the link checker, so a redirect onto them is ignored.
linkcheck_ignore = [
    # ACS returns 403 to the link checker.
    r"https://pubs\.acs\.org/.*",
    # Wiley returns 403 to the link checker.
    r"https://onlinelibrary\.wiley\.com/.*",
    # AIP returns 403 to the link checker.
    r"https://pubs\.aip\.org/.*",
    # Taylor & Francis returns 403 to the link checker.
    r"https://www\.tandfonline\.com/.*",
    # The dataset pages require a Hugging Face login.
    r"https://huggingface\.co/datasets/d-SEAMS/.*",
]
html_static_path = ["_static"]
html_css_files = ["custom.css"]
html_js_files = []
html_favicon = "_static/dSeamsLogo.ico"

# Mermaid: default CDN; diagrams via ``.. mermaid::`` from Org RST export.
mermaid_version = "11.4.0"
mermaid_height = "auto"
mermaid_init_js = (
    'import mermaid from "https://cdn.jsdelivr.net/npm/mermaid@11.4.0/dist/mermaid.esm.min.mjs";'
    "mermaid.initialize({startOnLoad:false, theme:'base', "
    "flowchart:{htmlLabels:true, padding:16, nodeSpacing:30, rankSpacing:36}, "
    "themeVariables:{fontSize:'18px', fontFamily:'Inter, sans-serif', "
    "primaryColor:'#b7e4d8', primaryTextColor:'#102220', primaryBorderColor:'#0c5c56', "
    "lineColor:'#0c5c56', secondaryColor:'#d5efe8', tertiaryColor:'#ffffff', "
    "background:'#ffffff'}});"
)

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
            "title": "C++ API",
            "url": "/api-cpp/index.html",
            "resource": True,
        },
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
                {
                    "title": "linkcell",
                    "url": "https://github.com/d-SEAMS/linkcell",
                    "summary": "Periodic linked-cell k-nearest neighbours",
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


def _copy_caption_tracks(app, exception):
    """Publish the draft VTT files next to the HTML so <track> can load them."""
    if exception is not None:
        return
    import shutil

    source = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "videos"))
    dest = os.path.join(app.outdir, "_static", "captions")
    os.makedirs(dest, exist_ok=True)
    if not os.path.isdir(source):
        return
    for name in os.listdir(source):
        track = os.path.join(source, name, "captions.vtt")
        if os.path.isfile(track):
            shutil.copy(track, os.path.join(dest, name + ".vtt"))


def setup(app):
    app.connect("build-finished", _copy_caption_tracks)
