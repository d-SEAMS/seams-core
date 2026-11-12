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
    ]
    OUT.mkdir(parents=True, exist_ok=True)
    parts = ["<html><body>"]
    for path in required:
        if not path.is_file():
            print("missing", path, file=sys.stderr)
            return 1
        body = path.read_text()
        html = ["<h1>{}</h1>".format(path.stem), "<pre>"]
        html.append(body.replace("&", "&amp;").replace("<", "&lt;"))
        html.append("</pre>")
        parts.extend(html)
    parts.append("</body></html>")
    (OUT / "index.html").write_text("\n".join(parts))
    print("docs ok", len(names), "symbols")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
