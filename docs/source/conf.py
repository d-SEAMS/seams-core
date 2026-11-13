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


def _write_binding_overhead():
    """Draw the binding-overhead figure before Sphinx reads the page."""
    import importlib.util

    script = os.path.join(_docs, "render_binding_overhead.py")
    spec = importlib.util.spec_from_file_location("render_binding_overhead", script)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    if module.main([]) != 0:
        raise SystemExit("binding overhead page failed")


_write_binding_overhead()


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


def _stage_markdown_book():
    """Copy the language pages into the Sphinx tree so myst renders them."""
    import pathlib

    groups = (
        ("api", "stable", "Stable API"),
        ("engines", "engines", "Engines"),
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
        body = [title, underline, "", ".. toctree::", "   :maxdepth: 1", ""]
        body.extend("   {}".format(name) for name in names)
        body.append("")
        (dest / "index.rst").write_text("\n".join(body), encoding="utf-8")


def _stage_ctypes_autodoc():
    """Pull numpydoc into the book.

    rgpot keeps the narrative page and the generated API apart. This page
    is the generated half. A missing import still leaves a page so the
    how-to link resolves.
    """
    stable = os.path.join(_source, "stable")
    page = os.path.join(stable, "seams_ctypes.rst")
    index = os.path.join(stable, "index.rst")
    try:
        import seams  # noqa: F401
    except Exception:
        body = (
            "ctypes reference\n"
            "================\n"
            "\n"
            "The ``seams`` package did not import in this build.\n"
            "The how-to is :doc:`python`.\n"
        )
    else:
        body = (
            "ctypes reference\n"
            "================\n"
            "\n"
            "Generated from the numpydoc strings on ``seams``.\n"
            "The how-to is :doc:`python`.\n"
            "\n"
            ".. automodule:: seams\n"
            "   :members:\n"
            "   :undoc-members:\n"
            "   :show-inheritance:\n"
        )
    with open(page, "w", encoding="utf-8") as handle:
        handle.write(body)
    if os.path.isfile(index):
        with open(index, "a", encoding="utf-8") as handle:
            handle.write("   seams_ctypes\n")


def _native_page(folder, title, detail):
    """One HTML page so a language link resolves when the tool is absent."""
    os.makedirs(folder, exist_ok=True)
    index = os.path.join(folder, "index.html")
    if os.path.isfile(index):
        return
    with open(index, "w", encoding="utf-8") as handle:
        handle.write("<!DOCTYPE html><title>{}</title><p>{}</p>\n".format(title, detail))


def _stage_native_html():
    """Copy Documenter and pkgdown HTML next to the Sphinx book.

    The language pages link at these paths. A missing tool writes one
    page and does not stop the book, the same way a missing C++ index
    still leaves the Reference toctree resolvable.
    """
    import glob
    import shutil
    import subprocess

    native = os.path.join(_docs, "native-html")
    os.makedirs(native, exist_ok=True)
    root = os.path.abspath(os.path.join(_docs, ".."))
    julia = shutil.which("julia")
    if not julia:
        matches = sorted(glob.glob(os.path.join(root, "bbdir", "julia-*", "bin", "julia")))
        julia = matches[-1] if matches else ""
    pairs = (
        ("julia-dseams", os.path.join(root, "bindings", "julia", "Dseams", "docs")),
        ("julia-seams", os.path.join(root, "bindings", "julia", "Seams", "docs")),
    )
    if julia:
        env = os.environ.copy()
        for name, project in pairs:
            subprocess.run(
                [julia, "--project=" + project, "-e", "using Pkg; Pkg.instantiate()"],
                cwd=root,
                env=env,
                check=False,
            )
            subprocess.run(
                [julia, "--project=" + project, os.path.join(project, "make.jl")],
                cwd=root,
                env=env,
                check=False,
            )
            built = os.path.join(project, "build")
            dest = os.path.join(native, name)
            if os.path.isdir(built):
                if os.path.isdir(dest):
                    shutil.rmtree(dest)
                shutil.copytree(built, dest)
                shutil.rmtree(built, ignore_errors=True)
            else:
                _native_page(dest, name, "Documenter did not write this build.")
    else:
        print("julia is absent; Documenter pages are placeholders")
        for name, _project in pairs:
            _native_page(os.path.join(native, name), name, "Julia is not on PATH.")
    r_dest = os.path.join(native, "pkgdown")
    rscript = shutil.which("Rscript")
    if rscript:
        env = os.environ.copy()
        libs = env.get("R_LIBS_USER", os.path.join(root, "bbdir", "rlib"))
        env["R_LIBS_USER"] = libs
        env.setdefault("DSEAMS_ROOT", root)
        done = subprocess.run(
            [
                rscript,
                "-e",
                "pkgdown::build_site('bindings/r/dseams', override=list(destination='%s'))" % r_dest,
            ],
            cwd=root,
            env=env,
            check=False,
        )
        if done.returncode != 0 or not os.path.isfile(os.path.join(r_dest, "index.html")):
            print("pkgdown did not build; the roxygen page is a placeholder")
            _native_page(r_dest, "dseams", "pkgdown did not build. The roxygen topics are in man/.")
    else:
        print("Rscript is absent; the pkgdown page is a placeholder")
        _native_page(r_dest, "dseams", "Rscript is not on PATH. The roxygen topics are in man/.")


def _link_capi():
    """Point the generated book at the committed C API section."""
    index = os.path.join(_source, "index.rst")
    if not os.path.isfile(index):
        return
    text = open(index, encoding="utf-8").read()
    if "capi/index" in text:
        return
    extra = (
        "\n.. toctree::\n"
        "   :maxdepth: 2\n"
        "   :caption: C API\n"
        "\n"
        "   capi/index\n"
    )
    with open(index, "a", encoding="utf-8") as handle:
        handle.write(extra)


_stage_markdown_book()
_stage_ctypes_autodoc()
_stage_native_html()
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
        "   :maxdepth: 2\n"
        "   :caption: Language APIs\n"
        "\n"
        "   lang/index\n"
    )
    with open(index, "a", encoding="utf-8") as handle:
        handle.write(extra)


_link_lang()
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

# Rust API pages. sphinxcontrib-rust (the same extension rgpot uses) renders
# rustdoc into the book when the safe crate can link libseams. Without that
# library, or without the extension, a short page keeps the toctree valid.
_rust_page = os.path.join(_source, "crates", "seams", "lib.rst")
_seams_lib_dir = os.environ.get("SEAMS_LIB_DIR", "")
_rust_ready = False
if _seams_lib_dir:
    try:
        import importlib

        importlib.import_module("sphinxcontrib_rust")
        _rust_ready = True
    except ImportError:
        _rust_ready = False
if _rust_ready:
    extensions.append("sphinxcontrib_rust")
    rust_crates = {
        "seams": os.path.abspath(os.path.join(_source, "../../bindings/rust/seams")),
    }
    rust_doc_dir = os.path.join(_source, "crates")
    rust_rustdoc_fmt = "rst"
    rust_generate_mode = "always"
else:
    os.makedirs(os.path.dirname(_rust_page), exist_ok=True)
    if not os.path.isfile(_rust_page):
        with open(_rust_page, "w", encoding="utf-8") as _rust_stub:
            _rust_stub.write(
                "Rust API\n"
                "========\n"
                "\n"
                "The safe crate is ``seams``. ``seams-sys`` is the raw FFI crate.\n"
                "``cargo doc --manifest-path bindings/rust/Cargo.toml --no-deps``\n"
                "writes the rustdoc. Set ``SEAMS_LIB_DIR`` and install\n"
                "``sphinxcontrib-rust`` to render those pages in this book.\n"
            )

templates_path = ["_templates"]
exclude_patterns = []
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
html_extra_path = [os.path.join(_docs, "native-html")]
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
            "title": "C++ API",
            "url": "/api-cpp/index.html",
            "resource": True,
        },
        {
            "title": "Fortran API",
            "url": "/fortran/index.html",
            "resource": True,
        },
        {
            "title": "Lua API",
            "url": "/lua/index.html",
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


def _skip_attributes(app, what, name, obj, skip, options):
    """Napoleon already documents dataclass fields from the Attributes section."""
    if what == "attribute":
        return True
    return skip


def setup(app):
    app.connect("build-finished", _copy_caption_tracks)
    app.connect("autodoc-skip-member", _skip_attributes)
