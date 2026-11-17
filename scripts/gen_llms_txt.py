#!/usr/bin/env python3
# Copyright (c) d-SEAMS developers
"""Write llms.txt, the plain-text map of the published d-SEAMS manual.

The pages and their order come from the toctrees of docs/orgmode/index.org,
one section per top-level caption. Nested toctrees in those pages are
walked in order, so a section index contributes its children. A page is
an Org file under docs/orgmode, a Markdown page under docs/api,
docs/engines or docs/tutorials (published as stable/, engines/ and
fronts/), or an rst page under docs/source. Each line links the published
HTML on https://docs.dseams.info with the page title and the first
sentence of its prose. A last section links the language binding
repositories.

llms.txt is committed at the repository root and Sphinx copies it to the
site root. ``--check`` exits 1 when the committed file differs from a
fresh generation.

>>> plain("Run ~seams rings~ on [[file:x.org][a frame]].")
'Run seams rings on a frame.'
>>> first_sentence("CHILL+ labels each oxygen by its bonds. It reads four.")
'CHILL+ labels each oxygen by its bonds.'
"""

from __future__ import annotations

import argparse
import posixpath
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DOCS = ROOT / "docs"
ORG = DOCS / "orgmode"
SITE = "https://docs.dseams.info"
OUTPUT = ROOT / "llms.txt"

TITLE = "d-SEAMS"
SUMMARY = (
    "d-SEAMS (seams-core) reads molecular-dynamics trajectories and names "
    "the ice in them: CHILL and CHILL+, topological cages, primitive rings, "
    "Steinhardt q_l, and ice and brine fronts, computed on one shared "
    "neighbour list. This repository builds the C++ engine, the C API "
    "libseams (include/seams.h), the C ABI libseams_lang that the language "
    "bindings link, the seams command, and the LAMMPS, GROMACS, OpenMM, "
    "PLUMED and in-situ engine drivers. The Python, Lua and Fennel, R, "
    "Julia, Rust and Fortran bindings live in their own repositories and "
    "release against the seams-core tag with the same version."
)

BINDINGS = (
    ("PydSEAMSlib", "https://github.com/d-SEAMS/PydSEAMSlib",
     "Python binding, pydseamslib on PyPI (import dseams, import pydseams)."),
    ("yodaStruct", "https://github.com/d-SEAMS/yodaStruct",
     "Lua and Fennel binding, the LuaRocks package dseams (require(\"dseams\"))."),
    ("dseams-r", "https://github.com/d-SEAMS/dseams-r",
     "R binding, the R package dseams."),
    ("Dseams.jl", "https://github.com/d-SEAMS/Dseams.jl",
     "Julia binding, the package Dseams over libseams_jll."),
    ("dseams-rs", "https://github.com/d-SEAMS/dseams-rs",
     "Rust binding, the crates dseams-sys and dseams."),
    ("dseams-fortran", "https://github.com/d-SEAMS/dseams-fortran",
     "Fortran binding, the fpm package dseams."),
)

COPIED = {
    "stable": DOCS / "api",
    "engines": DOCS / "engines",
    "fronts": DOCS / "tutorials",
}


def toctrees(text):
    """(caption, [docname, ...]) for each toctree in an Org index."""
    trees = []
    current = None
    for line in text.splitlines():
        stripped = line.strip()
        if stripped.startswith(".. toctree::"):
            current = ["", []]
            trees.append(current)
            continue
        if current is None:
            continue
        if stripped.startswith(":caption:"):
            current[0] = stripped.split(":", 2)[2].strip()
            continue
        if stripped.startswith(":") or stripped == "":
            continue
        if line.startswith("   ") and " " not in stripped and not stripped.startswith(".. "):
            current[1].append(stripped)
            continue
        current = None
    return [(caption, names) for caption, names in trees if names]


def _tracked_source_rst():
    """Rst pages that live in git. Sphinx writes the rest at build time."""
    listed = set()
    import subprocess

    result = subprocess.run(
        ["git", "ls-files", "-z", "--", "docs/source"],
        cwd=ROOT,
        check=True,
        capture_output=True,
    )
    for raw in result.stdout.split(b"\0"):
        if raw.endswith(b".rst"):
            listed.add(raw.decode())
    return listed


_TRACKED_RST = None


def source_of(name):
    global _TRACKED_RST
    org = ORG / (name + ".org")
    if org.is_file():
        return org
    head, _, rest = name.partition("/")
    if head in COPIED and rest:
        page = COPIED[head] / (rest + ".md")
        if page.is_file():
            return page
    rst = DOCS / "source" / (name + ".rst")
    if _TRACKED_RST is None:
        _TRACKED_RST = _tracked_source_rst()
    key = rst.relative_to(ROOT).as_posix()
    if key in _TRACKED_RST and rst.is_file():
        return rst
    return None


def plain(text):
    text = re.sub(r"\[\[(?:file:)?[^]\n]+\]\[([^]\n]+)\]\]", r"\1", text)
    text = re.sub(r"\[\[([^]\n]+)\]\]", r"\1", text)
    text = re.sub(r"`([^`<\n]+?)\s*<[^>\n]+>`_+", r"\1", text)
    text = re.sub(r"\[([^\]\n]+)\]\([^)\n]+\)", r"\1", text)
    text = re.sub(r":(?:doc|ref|math|cite[pt]?):`([^`<]*?)\s*(?:<[^>]*>)?`", r"\1", text)
    text = re.sub(r"~([^~\n]+)~", r"\1", text)
    text = re.sub(r"(?<![\w=])=([^=\n]+)=(?![\w=])", r"\1", text)
    text = text.replace("``", "").replace("`", "").replace("**", "")
    text = re.sub(r"\\[()\[\]]", "", text)
    return " ".join(text.split())


def first_sentence(text, limit=200):
    """Leading sentences of ``text`` until 30 characters, cut at ``limit``.

    A period after a short capitalised word (J. Phys. Chem.) is not an end.

    >>> first_sentence("Nguyen and Molinero, J. Phys. Chem. B 119, wrote the rule. It is old.")
    'Nguyen and Molinero, J. Phys. Chem. B 119, wrote the rule.'
    """
    taken = ""
    start = 0
    for match in re.finditer(r"[.!?](?=\s+|$)", text):
        word = text[: match.start()].rsplit(" ", 1)[-1]
        if re.fullmatch(r"[A-Z][a-z]{0,6}", word):
            continue
        taken = text[: match.end()].strip()
        start = match.end()
        if len(taken) >= 30:
            break
    text = taken if taken and start else text
    if len(text) > limit:
        text = text[: limit - 3].rsplit(" ", 1)[0] + "..."
    return text


def lead_paragraph(paragraphs):
    """The first paragraph with a sentence in it, not a bare version line."""
    for lines in paragraphs:
        text = re.sub(r"^d-SEAMS [0-9]+(?:\.[0-9]+)+\.\s*", "", plain(" ".join(lines)))
        if len(text) >= 25:
            return text
    return plain(" ".join(paragraphs[0])) if paragraphs else ""


def title_and_paragraph(path):
    text = path.read_text(encoding="utf-8")
    lines = text.splitlines()
    title = ""
    paragraphs = []
    paragraph = []
    skipping = False
    fence = False

    def flush():
        if paragraph:
            paragraphs.append(list(paragraph))
            paragraph.clear()

    if path.suffix == ".org":
        for line in lines:
            stripped = line.strip()
            upper = stripped.upper()
            if upper.startswith("#+TITLE:"):
                title = plain(stripped.split(":", 1)[1])
                continue
            if upper.startswith("#+BEGIN_"):
                flush()
                skipping = True
                continue
            if upper.startswith("#+END_"):
                skipping = False
                continue
            if skipping or not stripped or stripped.startswith(("#", "|", "*", "- ", ":")):
                flush()
                continue
            paragraph.append(stripped)
    elif path.suffix == ".md":
        for line in lines:
            stripped = line.strip()
            if stripped.startswith("```"):
                flush()
                fence = not fence
                continue
            if fence:
                continue
            if stripped.startswith("# ") and not title:
                title = plain(stripped[2:])
                continue
            if not stripped or stripped.startswith(("#", "<!--", "|", "-", ":::")):
                flush()
                continue
            paragraph.append(stripped)
    else:
        for index, line in enumerate(lines):
            stripped = line.strip()
            under = lines[index + 1].strip() if index + 1 < len(lines) else ""
            if not title and stripped and under and set(under) <= set("=-~^*#") and len(under) >= len(stripped):
                title = plain(stripped)
                continue
            if stripped and set(stripped) <= set("=-~^*#"):
                continue
            if stripped.startswith(".."):
                flush()
                skipping = True
                continue
            if skipping and (line.startswith(" ") or not stripped):
                continue
            skipping = False
            if not stripped or stripped.startswith(("*", "-", "|")):
                flush()
                continue
            if title:
                paragraph.append(stripped)
    flush()
    return title or path.stem, first_sentence(lead_paragraph(paragraphs))


def walk_pages(name, seen, items):
    """Append ``name`` and the pages named by its toctrees.

    Entries are relative to the page that lists them. ``..`` climbs out
    of that page's directory, the same way Sphinx resolves a toctree.
    """
    name = posixpath.normpath(name)
    if name in seen:
        return
    path = source_of(name)
    if path is None:
        return
    seen.add(name)
    title, summary = title_and_paragraph(path)
    entry = "- [{}]({}/{}.html)".format(title, SITE, name)
    items.append(entry + (": " + summary if summary else ""))
    if path.suffix != ".org":
        return
    parent = posixpath.dirname(name)
    for _caption, children in toctrees(path.read_text(encoding="utf-8")):
        for child in children:
            walk_pages(posixpath.join(parent, child), seen, items)


def render():
    index = (ORG / "index.org").read_text(encoding="utf-8")
    lines = ["# " + TITLE, "", "> " + SUMMARY, ""]
    lines.append(
        "Each link below is a page of the published manual on {}. "
        "Contributor instructions are on the Contributing and Repository "
        "map pages.".format(SITE)
    )
    lines.append("")
    seen = set()
    home_title, home_summary = title_and_paragraph(ORG / "index.org")
    lines.extend(["## Start", "", "- [{}]({}/index.html): {}".format(home_title, SITE, home_summary), ""])
    seen.add("index")
    for caption, names in toctrees(index):
        items = []
        for name in names:
            walk_pages(name, seen, items)
        if items:
            lines.extend(["## " + caption, ""])
            lines.extend(items)
            lines.append("")
    lines.extend(["## Language bindings", ""])
    for name, url, summary in BINDINGS:
        lines.append("- [{}]({}): {}".format(name, url, summary))
    return "\n".join(lines).rstrip() + "\n"


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--check", action="store_true", help="exit 1 when llms.txt differs from a fresh generation")
    args = parser.parse_args(argv)
    text = render()
    if args.check:
        current = OUTPUT.read_text(encoding="utf-8") if OUTPUT.is_file() else None
        if current != text:
            print("llms.txt is out of date. Run: python3 scripts/gen_llms_txt.py", file=sys.stderr)
            return 1
        print("llms.txt ok ({} links)".format(text.count("\n- [")))
        return 0
    OUTPUT.write_text(text, encoding="utf-8")
    print("wrote", OUTPUT.relative_to(ROOT))
    return 0


if __name__ == "__main__":
    sys.exit(main())
