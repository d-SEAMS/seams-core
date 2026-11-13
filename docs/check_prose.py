#!/usr/bin/env python3
# Copyright (c) d-SEAMS developers
"""Fail when public docs use the filler this release does not.

The scan is docs/**/*.md, docs/**/*.rst, README.md, and changelog.d.
Command flags and fenced examples are left alone. A double hyphen with
spaces on both sides, outside those examples, is punctuation.
"""

import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DOCS = ROOT / "docs"
SKIP = {"build", "doctrees", "__pycache__"}
WORDS = re.compile(
    r"\b(delve|delves|delving|leverage|leverages|leveraging|leveraged|"
    r"robust|seamless|underpin|underpins|underpinning|"
    r"end-to-end|gate|gates|shipping)\b",
    re.I,
)
TRUTH = re.compile(r"the truth is", re.I)
NOT_BUT = re.compile(r"\bnot\b[^.!?\n]{0,80}\bbut\b", re.I)
EM = "\u2014"
PUNCT = re.compile(r"\s--\s")
NIX = re.compile(r"\bnix run\b")
CARGO_SPLIT = re.compile(r" -- --")


def pages():
    found = []
    for path in DOCS.rglob("*"):
        if not path.is_file() or path.suffix not in {".md", ".rst"}:
            continue
        if SKIP.intersection(path.relative_to(DOCS).parts):
            continue
        found.append(path)
    readme = ROOT / "README.md"
    if readme.is_file():
        found.append(readme)
    found.extend(sorted((ROOT / "changelog.d").glob("*.md")))
    return found


def in_fence(line, state):
    """Track markdown fences, Org source blocks, and reST code blocks."""
    stripped = line.strip()
    if stripped.startswith("```") or stripped.startswith("~~~"):
        state["fence"] = not state["fence"]
        return True
    if stripped.startswith("#+begin_src") or stripped.startswith("#+BEGIN_SRC"):
        state["org"] = True
        return True
    if stripped.startswith("#+end_src") or stripped.startswith("#+END_SRC"):
        state["org"] = False
        return True
    if state["org"] or state["fence"]:
        return True
    if stripped.startswith(".. code") or stripped.endswith("::"):
        state["rest"] = True
        state["rest_indent"] = None
        return True
    if state["rest"]:
        if stripped == "":
            return True
        indent = len(line) - len(line.lstrip(" "))
        if state["rest_indent"] is None:
            if indent >= 3:
                state["rest_indent"] = indent
                return True
            state["rest"] = False
            return False
        if indent >= state["rest_indent"]:
            return True
        state["rest"] = False
    return False


def punctuation(line):
    if not PUNCT.search(line):
        return False
    if NIX.search(line) or CARGO_SPLIT.search(line):
        return False
    return True


def main():
    errors = []
    for path in pages():
        state = {"fence": False, "org": False, "rest": False, "rest_indent": None}
        for number, line in enumerate(path.read_text(errors="replace").splitlines(), start=1):
            if in_fence(line, state):
                continue
            rel = path.relative_to(ROOT)
            if WORDS.search(line):
                errors.append("{}:{} {}".format(rel, number, line.strip()[:120]))
            elif TRUTH.search(line):
                errors.append("{}:{} {}".format(rel, number, line.strip()[:120]))
            elif NOT_BUT.search(line):
                errors.append("{}:{} {}".format(rel, number, line.strip()[:120]))
            elif EM in line:
                errors.append("{}:{} {}".format(rel, number, line.strip()[:120]))
            elif punctuation(line):
                errors.append("{}:{} {}".format(rel, number, line.strip()[:120]))
    if errors:
        print("\n".join(errors), file=sys.stderr)
        print("{} prose hits".format(len(errors)), file=sys.stderr)
        return 1
    print("prose check ok")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
