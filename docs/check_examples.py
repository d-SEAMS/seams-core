#!/usr/bin/env python3
# Copyright (c) d-SEAMS developers
"""Fail when a documented example is not a slice of a file the tests run.

Markdown pages under docs/api, docs/engines, and docs/tutorials are copied
into docs/source before Sphinx runs. A literalinclude path is resolved from
that copy. An Org or Markdown fence whose first line is

    # excerpt: path/from/repo/root

must be followed by text that occurs, unchanged, in that file.
"""

import pathlib
import re
import sys

ROOT = pathlib.Path(__file__).resolve().parents[1]
DOCS = ROOT / "docs"

COPY_ROOT = {
    "api": DOCS / "source" / "stable",
    "engines": DOCS / "source" / "engines",
    "tutorials": DOCS / "source" / "fronts",
}

DIRECTIVE = re.compile(
    r"^(?:#|//|!)\s*excerpt:\s*(\S+)\s*$"
)
ORG_SRC = re.compile(r"#\+begin_src[^\n]*\n(.*?)#\+end_src", re.S | re.I)


def c_symbols():
    text = (ROOT / "include" / "seams.h").read_text()
    found = re.findall(r"\bseams_[a-z0-9_]+\b", text)
    skip = {
        "seams_h",
        "seams_system",
        "seams_neighbours",
        "seams_rings",
        "seams_status",
        "seams_label",
        "seams_ring_symmetry",
        "seams_chill_rule",
        "seams_graph",
        "seams_version",
    }
    return sorted({name for name in found if name not in skip and not name.endswith("_t")})


def lang_symbols():
    text = (ROOT / "bindings" / "cxx" / "seams_lang.h").read_text()
    found = re.findall(r"\b(seams_[a-z0-9_]+)\s*\(", text)
    seen = []
    for name in found:
        if name not in seen:
            seen.append(name)
    return seen


def fortran_symbols():
    text = (ROOT / "bindings" / "fortran" / "seams.f90").read_text()
    return re.findall(r'bind\(C, name="([^"]+)"\)', text)


def present(text, token):
    return re.search(r"(?<!\w)" + re.escape(token) + r"(?!\w)", text) is not None


def parse_lines(spec, nlines):
    chosen = []
    for part in spec.split(","):
        part = part.strip()
        if not part:
            continue
        if "-" in part:
            start, end = part.split("-", 1)
            lo, hi = int(start), int(end)
        else:
            lo = hi = int(part)
        if lo < 1 or hi > nlines or hi < lo:
            raise ValueError("bad line range {} (file has {})".format(part, nlines))
        chosen.extend(range(lo, hi + 1))
    return chosen


def copied_dir(path):
    rel = path.relative_to(DOCS)
    return COPY_ROOT[rel.parts[0]]


def resolve_include(doc, raw):
    """raw is the path Sphinx sees from the copied page."""
    base = copied_dir(doc) / doc.name
    return (base.parent / raw).resolve()


def load_lines(path):
    return path.read_text().splitlines()


class Corpus:
    def __init__(self):
        self.by_doc = {}

    def add(self, doc, text):
        self.by_doc.setdefault(doc, []).append(text)

    def text(self, doc):
        return "\n".join(self.by_doc.get(doc, []))


def harvest_markdown(doc, raw, errors, corpus):
    lines = raw.splitlines()
    index = 0
    while index < len(lines):
        line = lines[index]
        if line.startswith("```{literalinclude}"):
            parts = line.split()
            if len(parts) < 2:
                errors.append("{}: literalinclude has no path".format(doc))
                index += 1
                continue
            include = parts[1]
            index += 1
            options = []
            while index < len(lines) and not lines[index].startswith("```"):
                options.append(lines[index])
                index += 1
            if index < len(lines) and lines[index].startswith("```"):
                index += 1
            line_opt = None
            for option in options:
                stripped = option.strip()
                if stripped.startswith(":lines:"):
                    line_opt = stripped.split(":", 2)[-1].strip()
            try:
                target = resolve_include(doc, include)
            except Exception as exc:  # noqa: BLE001
                errors.append("{}: {}".format(doc, exc))
                continue
            if not target.is_file():
                errors.append("{}: missing {}".format(doc, target))
                continue
            source_lines = load_lines(target)
            if line_opt:
                try:
                    chosen = parse_lines(line_opt, len(source_lines))
                except ValueError as exc:
                    errors.append("{}: {}".format(doc, exc))
                    continue
                snippet = "\n".join(source_lines[i - 1] for i in chosen)
            else:
                snippet = "\n".join(source_lines)
            corpus.add(doc, snippet)
            continue
        if line.startswith("```"):
            index += 1
            body = []
            while index < len(lines) and not lines[index].startswith("```"):
                body.append(lines[index])
                index += 1
            if index < len(lines) and lines[index].startswith("```"):
                index += 1
            record_fence(doc, "\n".join(body), errors, corpus)
            continue
        index += 1


def harvest(errors):
    corpus = Corpus()
    pages = []
    for folder in ("api", "engines", "tutorials"):
        pages.extend(sorted((DOCS / folder).glob("*.md")))
    pages.extend(sorted((DOCS / "orgmode").rglob("*.org")))
    for doc in pages:
        raw = doc.read_text()
        if doc.suffix == ".md":
            harvest_markdown(doc, raw, errors, corpus)
        else:
            for match in ORG_SRC.finditer(raw):
                record_fence(doc, match.group(1), errors, corpus)
    return corpus


def record_fence(doc, body, errors, corpus):
    lines = body.splitlines()
    if not lines:
        return
    directive = DIRECTIVE.match(lines[0].strip())
    if directive is None:
        return
    target = (ROOT / directive.group(1)).resolve()
    excerpt = "\n".join(lines[1:]).strip("\n")
    if not target.is_file():
        errors.append("{}: excerpt file missing {}".format(doc, target))
        return
    source = target.read_text()
    if excerpt not in source:
        errors.append("{}: excerpt is not in {}".format(doc, directive.group(1)))
        return
    corpus.add(doc, excerpt)


def require_tokens(corpus, docs, tokens, errors, label):
    blob = "\n".join(corpus.text(doc) for doc in docs)
    missing = [token for token in tokens if not present(blob, token)]
    if missing:
        names = ", ".join(str(doc.relative_to(ROOT)) for doc in docs)
        errors.append("{} missing from {}: {}".format(label, names, " ".join(missing)))


def main():
    errors = []
    corpus = harvest(errors)
    api = DOCS / "api"
    engines = DOCS / "engines"
    howto = DOCS / "orgmode" / "howto"
    require_tokens(corpus, [api / "c.md"], c_symbols(), errors, "C API")
    require_tokens(corpus, [api / "cpp.md"], lang_symbols(), errors, "language ABI")
    require_tokens(corpus, [api / "fortran.md"], fortran_symbols(), errors, "Fortran")
    require_tokens(
        corpus,
        [api / "python.md"],
        [
            "version",
            "set_num_threads",
            "num_threads",
            "Frame",
            "set_frame",
            "count",
            "chill_plus",
            "q4q6",
            "neighbours",
            "rings",
            "cages",
            "signature_cages",
            "symmetry",
            "from_ase",
            "from_mdanalysis",
            "PUBLISHED",
            "OVITO",
            "CUTOFF",
        ],
        errors,
        "Python",
    )
    require_tokens(
        corpus,
        [api / "julia.md"],
        [
            "read_frame",
            "nop",
            "nframes",
            "chill_plus",
            "cages",
            "primitive_rings",
            "signature_cages",
            "symmetry",
            "steinhardt",
            "q4_q6",
            "jump_tau90",
            "set_num_threads",
            "num_threads",
            "version",
            "set_ortho!",
            "q4q6",
            "six_rings",
            "System",
        ],
        errors,
        "Julia",
    )
    require_tokens(
        corpus,
        [api / "r.md"],
        [
            "read_frame",
            "positions",
            "chill_plus",
            "cages",
            "primitive_rings",
            "signature_cages",
            "symmetry",
            "steinhardt",
            "nframes",
            "print",
            "q4_q6",
            "jump_tau90",
        ],
        errors,
        "R",
    )
    require_tokens(
        corpus,
        [api / "rust.md"],
        [
            "set_num_threads",
            "num_threads",
            "version",
            "System::new",
            "set_ortho",
            "set_frame",
            "count",
            "chill_plus",
            "q4q6",
            "neighbour_degree",
            "six_ring_fraction",
        ],
        errors,
        "Rust",
    )
    require_tokens(
        corpus,
        [engines / "lammps.md"],
        ["fix ext all external", "compute ice all dseams", "seams_chill_plus"],
        errors,
        "LAMMPS",
    )
    require_tokens(
        corpus,
        [engines / "gromacs.md"],
        ["gromacs_box_to_lammps", "seams_chill_plus", "gmx mdrun"],
        errors,
        "GROMACS",
    )
    require_tokens(
        corpus,
        [engines / "openmm.md"],
        ["DseamsReporter", "seams_chill_plus"],
        errors,
        "OpenMM",
    )
    require_tokens(
        corpus,
        [engines / "plumed.md"],
        ["SEAMS", "seams_q4q6", "seams_rings_find"],
        errors,
        "PLUMED",
    )
    require_tokens(
        corpus,
        [engines / "insitu.md"],
        ["observe", "seams_frame_chill_rings"],
        errors,
        "in-situ",
    )
    require_tokens(
        corpus,
        [howto / "rings.org"],
        ["rings"],
        errors,
        "rings how-to",
    )
    require_tokens(
        corpus,
        [howto / "cages-trace.org"],
        ["cages", "signature"],
        errors,
        "cages how-to",
    )
    require_tokens(
        corpus,
        [howto / "q4-q6.org"],
        ["steinhardt"],
        errors,
        "q4/q6 how-to",
    )
    require_tokens(
        corpus,
        [howto / "chill-plus.org"],
        ["chill-plus"],
        errors,
        "CHILL+ how-to",
    )
    require_tokens(
        corpus,
        [howto / "symmetry-classifier.org"],
        ["seams_frame_symmetry"],
        errors,
        "symmetry how-to",
    )
    require_tokens(
        corpus,
        [howto / "linkcell-carried.org"],
        ["orthoFullCsr", "carried"],
        errors,
        "linkcell how-to",
    )
    if errors:
        print("\n".join(errors), file=sys.stderr)
        return 1
    print("docs examples ok")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
