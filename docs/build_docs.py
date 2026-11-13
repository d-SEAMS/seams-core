# Copyright (c) d-SEAMS developers
"""Render the API pages and check that every C symbol is named."""

import pathlib
import re
import sys

ROOT = pathlib.Path(__file__).resolve().parents[1]
HEADER = ROOT / "include" / "seams.h"
DOCS = ROOT / "docs"
OUT = DOCS / "_build"


def symbols():
    text = HEADER.read_text()
    found = re.findall(r"\bseams_[a-z0-9_]+\b", text)
    skip = {"seams_h", "seams_system", "seams_neighbours", "seams_rings", "seams_status",
            "seams_label", "seams_ring_symmetry", "seams_chill_rule", "seams_graph",
            "seams_version"}
    return sorted({name for name in found if name not in skip and not name.endswith("_t")})


def render_markdown(text):
    """HTML for headings, lists, and fenced code. Code fences stay preformatted."""
    out = []
    paragraph = []
    in_code = False
    list_open = False

    def inline(raw):
        escaped = raw.replace("&", "&amp;").replace("<", "&lt;")
        parts = escaped.split("`")
        body = []
        for index, part in enumerate(parts):
            if index % 2:
                body.append("<code>{}</code>".format(part))
            else:
                body.append(part)
        return "".join(body)

    def flush_paragraph():
        nonlocal paragraph
        if paragraph:
            out.append("<p>{}</p>".format(inline(" ".join(paragraph))))
            paragraph = []

    def close_list():
        nonlocal list_open
        if list_open:
            out.append("</ul>")
            list_open = False

    for line in text.splitlines():
        if line.startswith("```"):
            flush_paragraph()
            close_list()
            if in_code:
                out.append("</code></pre>")
                in_code = False
            else:
                out.append("<pre><code>")
                in_code = True
            continue
        if in_code:
            out.append(line.replace("&", "&amp;").replace("<", "&lt;"))
            continue
        if not line.strip():
            flush_paragraph()
            close_list()
            continue
        if line.startswith("#"):
            flush_paragraph()
            close_list()
            level = len(line) - len(line.lstrip("#"))
            level = min(max(level, 1), 4)
            out.append("<h{0}>{1}</h{0}>".format(level, inline(line[level:].strip())))
            continue
        if line.startswith(("- ", "* ")):
            flush_paragraph()
            if not list_open:
                out.append("<ul>")
                list_open = True
            out.append("<li>{}</li>".format(inline(line[2:].strip())))
            continue
        close_list()
        paragraph.append(line.strip())
    flush_paragraph()
    close_list()
    if in_code:
        out.append("</code></pre>")
    return "\n".join(out)


def main():
    names = symbols()
    page = (DOCS / "api" / "c.md").read_text()
    missing = [name for name in names if name not in page]
    if missing:
        print("missing from docs/api/c.md:", " ".join(missing), file=sys.stderr)
        return 1
    required = [
        DOCS / "quickstart.md",
        DOCS / "api" / "python.md",
        DOCS / "api" / "cpp.md",
        DOCS / "api" / "julia.md",
        DOCS / "api" / "rust.md",
        DOCS / "api" / "fortran.md",
        DOCS / "api" / "r.md",
        DOCS / "engines" / "lammps.md",
        DOCS / "engines" / "gromacs.md",
        DOCS / "engines" / "openmm.md",
        DOCS / "engines" / "plumed.md",
        DOCS / "engines" / "insitu.md",
        DOCS / "tutorials" / "c.md",
        DOCS / "tutorials" / "cpp.md",
        DOCS / "tutorials" / "python.md",
        DOCS / "tutorials" / "julia.md",
        DOCS / "tutorials" / "rust.md",
        DOCS / "tutorials" / "fortran.md",
        DOCS / "tutorials" / "r.md",
        DOCS / "tutorials" / "lammps.md",
        DOCS / "tutorials" / "gromacs.md",
        DOCS / "tutorials" / "openmm.md",
        DOCS / "tutorials" / "plumed.md",
        DOCS / "tutorials" / "insitu.md",
        DOCS / "tutorials" / "engines.md",
    ]
    OUT.mkdir(parents=True, exist_ok=True)
    parts = ["<html><body>"]
    for path in required:
        if not path.is_file():
            print("missing", path, file=sys.stderr)
            return 1
        body = path.read_text()
        parts.append(render_markdown(body))
    parts.append("</body></html>")
    (OUT / "index.html").write_text("\n".join(parts))
    import importlib.util

    sitting = DOCS / "render_sitting.py"
    spec = importlib.util.spec_from_file_location("render_sitting", sitting)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    code = module.main(["--check"])
    if code != 0:
        return code
    checker = DOCS / "check_examples.py"
    spec = importlib.util.spec_from_file_location("check_examples", checker)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    code = module.main()
    if code != 0:
        return code
    print("docs ok", len(names), "symbols")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
