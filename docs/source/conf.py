#!/usr/bin/env python3

import multiprocessing
import os
import sys

# The science gallery loads OVITO, which leaves a Qt thread pool. sphinx-sitemap
# then forks a manager. Qt's atfork handler deadlocks that child, so the
# manager has to be a fresh interpreter.
try:
    multiprocessing.set_start_method("spawn")
except RuntimeError:
    pass

# -- Path setup --------------------------------------------------------------
_source = os.path.dirname(os.path.abspath(__file__))
_docs = os.path.dirname(_source)


def _stage_records():
    """Copy the frozen benchmark tables and figures into the Sphinx tree.

    The d-SEAMS 2.0 benchmark records moved to a dataset; docs/records
    holds the tables rendered from them. docs/records/README.md names
    the dataset and revision.
    """
    import shutil

    records = os.path.join(_docs, "records")
    for folder, _dirs, files in os.walk(records):
        rel = os.path.relpath(folder, records)
        for name in files:
            if name == "README.md":
                continue
            dest = os.path.join(_source, rel, name)
            os.makedirs(os.path.dirname(dest), exist_ok=True)
            shutil.copyfile(os.path.join(folder, name), dest)


_stage_records()


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


def _docs_ci():
    return os.environ.get("DOCS_CI") == "1"


def _sphinx_builder():
    """The -b name. Sphinx defaults to html when the flag is absent."""
    argv = sys.argv
    if "-b" in argv:
        index = argv.index("-b")
        if index + 1 < len(argv):
            return argv[index + 1]
    return "html"


def _publish_build():
    """DOCS_CI fails the book when a generated page is missing.

    linkcheck and doctest read the sources that job builds. They do not
    render molecule stills or the native API HTML, so a missing picture
    stays a publish error and is not an error for those builders.
    """
    return _docs_ci() and _sphinx_builder() not in ("linkcheck", "doctest")


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
    if status in (0, None):
        return
    if _sphinx_builder() in ("linkcheck", "doctest"):
        print("sphinx: science gallery is not part of this builder", file=sys.stderr)
        return
    kept = os.path.isdir(os.path.join(_source, "tutorials", "gallery"))
    if _docs_ci() or not kept:
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
    if status in (0, None):
        return
    if _sphinx_builder() in ("linkcheck", "doctest"):
        print("sphinx: molecule stills are not part of this builder", file=sys.stderr)
        return
    pictures = os.path.join(_source, "tutorials", "molecules")
    if _docs_ci() or not os.path.isdir(pictures):
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
    module.write_topic_embeds(videos, pathlib.Path(_source) / "videos" / "embeds")
    static = os.path.join(_source, "_static", "captions")
    os.makedirs(static, exist_ok=True)
    for name in os.listdir(videos):
        track = os.path.join(videos, name, "captions.vtt")
        if os.path.isfile(track):
            shutil.copy(track, os.path.join(static, name + ".vtt"))


_write_video_players()


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


def _ensure_root_index():
    """Give Sphinx a landing page when the Org export has not run yet.

    The book export writes this file. Linkcheck still walks every page
    the staging steps have already copied.
    """
    index = os.path.join(_source, "index.rst")
    if os.path.isfile(index):
        return
    sections = (
        "stable/index",
        "engines/index",
        "fronts/index",
        "capi/index",
        "lang/index",
        "science/index",
        "walkthrough/index",
        "api/index",
    )
    present = [name for name in sections if os.path.isfile(os.path.join(_source, name + ".rst"))]
    lines = ["d-SEAMS", "=======", "", ".. toctree::", "   :maxdepth: 2", ""]
    for name in present:
        lines.append("   " + name)
    lines.append("")
    with open(index, "w", encoding="utf-8") as handle:
        handle.write("\n".join(lines))


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


_stage_markdown_book()
_ensure_root_index()
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


def _toctree_docnames():
    """Docnames already named by a ``.. toctree::`` entry."""
    included = set()
    for dirpath, dirnames, filenames in os.walk(_source):
        dirnames[:] = [name for name in dirnames if name not in {"_static", "_templates"}]
        for filename in filenames:
            if not filename.endswith(".rst"):
                continue
            parent = os.path.relpath(dirpath, _source)
            parent = "" if parent == "." else parent.replace(os.sep, "/")
            lines = open(os.path.join(dirpath, filename), encoding="utf-8", errors="replace")
            in_tree = False
            for raw in lines:
                if raw.strip() == ".. toctree::":
                    in_tree = True
                    continue
                if not in_tree:
                    continue
                if raw.strip() == "":
                    continue
                if not raw.startswith("   "):
                    in_tree = False
                    continue
                entry = raw.strip()
                if entry.startswith(":"):
                    continue
                name = entry.split()[0].lstrip("/")
                if parent:
                    name = parent + "/" + name
                included.add(name)
    return included


def _include_orphan_pages():
    """Put every staged page in one toctree so linkcheck reads its links.

    Section indexes list some pages with ``:doc:`` only. Those pages are
    still part of the book, and a warning-as-error build would stop
    before the URL check.
    """
    index = os.path.join(_source, "index.rst")
    if not os.path.isfile(index):
        return
    start = ".. linkcheck-pages-start"
    end = ".. linkcheck-pages-end"
    text = open(index, encoding="utf-8").read()
    if start in text and end in text:
        text = text[: text.index(start)] + text[text.index(end) + len(end) :]
        text = text.rstrip() + "\n"
        with open(index, "w", encoding="utf-8") as handle:
            handle.write(text)
    included = _toctree_docnames()
    skip_prefixes = ("science/figures/", "xml/", "doxygen-include/")
    orphans = []
    for dirpath, dirnames, filenames in os.walk(_source):
        dirnames[:] = [name for name in dirnames if name not in {"_static", "_templates"}]
        for filename in filenames:
            if not filename.endswith((".rst", ".md")):
                continue
            rel = os.path.relpath(os.path.join(dirpath, filename), _source).replace(os.sep, "/")
            if rel == "index.rst" or rel.startswith(skip_prefixes):
                continue
            docname = os.path.splitext(rel)[0]
            if docname not in included:
                orphans.append(docname)
    if not orphans:
        return
    lines = [text.rstrip(), "", start, "", ".. toctree::", "   :hidden:", "   :maxdepth: 1", ""]
    for docname in sorted(orphans):
        lines.append("   " + docname)
    lines.extend(["", end, ""])
    with open(index, "w", encoding="utf-8") as handle:
        handle.write("\n".join(lines))


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


def _stub_unread_images():
    """linkcheck and doctest do not draw the molecule stills.

    The pages still name the pictures. A readable file keeps ``-W`` from
    failing the book on an image that this builder is not supposed to render.
    """
    if _sphinx_builder() not in ("linkcheck", "doctest"):
        return
    import pathlib
    import re

    source = pathlib.Path(_source).resolve()
    pattern = re.compile(r"(?<![\w/])((?:\.\./)*[A-Za-z0-9_./-]+\.(?:png|svg|jpg))")
    wanted = set()
    for page in source.rglob("*"):
        if page.suffix.lower() not in {".md", ".rst"}:
            continue
        try:
            body = page.read_text(encoding="utf-8", errors="replace")
        except OSError:
            continue
        for match in pattern.findall(body):
            target = (page.parent / match).resolve()
            try:
                target.relative_to(source)
            except ValueError:
                continue
            wanted.add(target)
    for target in sorted(wanted):
        if target.is_file():
            continue
        target.parent.mkdir(parents=True, exist_ok=True)
        if target.suffix.lower() == ".png":
            _write_png_stub(str(target))
        elif target.suffix.lower() == ".svg":
            target.write_text(
                '<svg xmlns="http://www.w3.org/2000/svg" width="1" height="1"/>\n',
                encoding="utf-8",
            )
        else:
            target.write_bytes(b"")


_stub_unread_images()
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

_include_orphan_pages()

# -- Project information -----------------------------------------------------
project = "d-SEAMS"
copyright = "2019-2026, d-SEAMS core team"
author = "Rohit Goswami, Amrita Goswami, Ruhila S"
version = "2.14.0"
release = "2.14.0"

# -- General configuration ---------------------------------------------------
extensions = [
    "myst_parser",
    "sphinx.ext.autodoc",
    "sphinx.ext.doctest",
    "sphinx.ext.intersphinx",
    "sphinx.ext.napoleon",
    "numpydoc",
    "sphinx_sitemap",
    "sphinx_design",
    "sphinxcontrib.mermaid",
]
napoleon_numpy_docstring = True
napoleon_google_docstring = False
napoleon_include_special_with_doc = True
autodoc_typehints = "description"
sys.path.insert(0, os.path.abspath(os.path.join(_docs, "..", "bindings", "python", "src")))
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
    # APS returns 403 to the link checker.
    r"https://journals\.aps\.org/.*",
    # RSC returns 403 to the link checker.
    r"https://pubs\.rsc\.org/.*",
    # ScienceDirect returns 403 to the link checker.
    r"https://www\.sciencedirect\.com/.*",
    # The dataset pages require a Hugging Face login.
    r"https://huggingface\.co/datasets/d-SEAMS/.*",
]
html_static_path = ["_static"]
# llms.txt at the site root: a plain-text map of the published pages.
html_extra_path = [os.path.join(os.path.dirname(_docs), "llms.txt")]
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
                    "title": "PydSEAMSlib (Python)",
                    "url": "https://d-seams.github.io/PydSEAMSlib/",
                    "summary": "pydseamslib: import dseams and import pydseams",
                    "external": True,
                },
                {
                    "title": "yodaStruct (Lua, Fennel)",
                    "url": "https://d-seams.github.io/yodaStruct/",
                    "summary": "require(\"dseams\")",
                    "external": True,
                },
                {
                    "title": "dseams-r (R)",
                    "url": "https://github.com/d-SEAMS/dseams-r",
                    "summary": "The R package dseams",
                    "external": True,
                },
                {
                    "title": "Dseams.jl (Julia)",
                    "url": "https://github.com/d-SEAMS/Dseams.jl",
                    "summary": "The Julia package Dseams",
                    "external": True,
                },
                {
                    "title": "dseams-rs (Rust)",
                    "url": "https://github.com/d-SEAMS/dseams-rs",
                    "summary": "The crates dseams-sys and dseams",
                    "external": True,
                },
                {
                    "title": "dseams-fortran (Fortran)",
                    "url": "https://github.com/d-SEAMS/dseams-fortran",
                    "summary": "The fpm package dseams",
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

# Shibuya writes the canonical link and og:url from this origin (pageurl).
# The trailing slash stays: sphinx-sitemap joins html_baseurl to the page
# link. Its default scheme is "{lang}{version}{link}". Sphinx sets language
# to "en", so that scheme emits https://docs.dseams.info/en/... and those
# paths 404. The HTML files live at the site root.
html_baseurl = "https://docs.dseams.info/"
sitemap_url_scheme = "{link}"


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


def _skip_attributes(app, what, name, obj, skip, options):
    """Napoleon already documents dataclass fields from the Attributes section."""
    if what == "attribute":
        return True
    return skip


def _relative_rst_links(app, doctree, docname):
    """Org export writes ``file.org`` as ``name.rst``. Publish those as HTML."""
    from docutils import nodes

    for node in doctree.findall(nodes.reference):
        uri = node.get("refuri")
        if not uri or "://" in uri or uri.startswith(("#", "mailto:")):
            continue
        path, sep, frag = uri.partition("#")
        if not path.endswith(".rst"):
            continue
        node["refuri"] = path[:-4] + ".html" + (sep + frag if sep else "")


def _edit_org_source(app, pagename, templatename, context, doctree):
    """Point Edit this page at the Org file on d-SEAMS/seams-core.

    Generated ``.rst`` is not in git. The repository is ``html_context``,
    the same owner and repo Shibuya uses for the sidebar link.
    """

    def edit_source_link(filename):
        if not isinstance(filename, str) or not filename.endswith(".rst"):
            return None
        root = os.path.abspath(os.path.join(app.srcdir, "..", ".."))
        rel = os.path.join("docs", "orgmode", filename[:-4] + ".org")
        if not os.path.isfile(os.path.join(root, rel)):
            return None
        return "https://github.com/{}/{}/blob/{}/{}".format(
            html_context["source_user"],
            html_context["source_repo"],
            html_context["source_version"],
            rel.replace(os.sep, "/"),
        )

    context["edit_source_link"] = edit_source_link


def setup(app):
    app.connect("doctree-resolved", _relative_rst_links)
    app.connect("html-page-context", _edit_org_source)
    app.connect("build-finished", _copy_caption_tracks)
    app.connect("autodoc-skip-member", _skip_attributes)
