#!/usr/bin/env python3
"""Inventory Lean declarations in lean/.

A declaration is a top-level ``theorem``, ``lemma``, ``def``, ``abbrev``,
``instance``, ``structure``, ``inductive`` or ``axiom``, optionally
prefixed by ``private``, ``protected`` or ``noncomputable``. The paper
claimed 55 declarations (53 public, 2 private) and 34 printed names.
Those totals are not stored here; this script recounts the sources.

    python3 scripts/count_lean.py [LEAN_DIR]
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

KINDS = ("theorem", "lemma", "def", "abbrev", "instance", "structure", "inductive", "axiom")
MODIFIERS = {"private", "protected", "noncomputable", "unsafe"}


def repo_root() -> Path:
    return Path(__file__).resolve().parents[1]


def declarations(path: Path) -> list[dict]:
    found = []
    for lineno, raw in enumerate(path.read_text(encoding="utf-8").splitlines(), 1):
        line = raw.strip()
        if not line or line.startswith("--") or line.startswith("/-") or line.startswith("-"):
            continue
        parts = line.split()
        private = False
        i = 0
        while i < len(parts) and parts[i] in MODIFIERS:
            if parts[i] == "private":
                private = True
            i += 1
        if i >= len(parts) or parts[i] not in KINDS:
            continue
        kind = parts[i]
        name = parts[i + 1] if i + 1 < len(parts) else ""
        name = name.split(":")[0].split("(")[0].split("{")[0].split("[")[0]
        if not name or name in {"where", ":="}:
            continue
        found.append(
            {
                "file": str(path),
                "line": lineno,
                "kind": kind,
                "name": name,
                "private": private,
            }
        )
    return found


def main() -> int:
    root = repo_root()
    lean = Path(sys.argv[1] if len(sys.argv) > 1 else root / "lean")
    if not lean.is_absolute():
        lean = (Path.cwd() / lean).resolve()
    files = sorted(lean.rglob("*.lean"))
    rows = []
    for path in files:
        rows.extend(declarations(path))
    public = [r for r in rows if not r["private"]]
    private = [r for r in rows if r["private"]]
    theorems = [r for r in rows if r["kind"] in {"theorem", "lemma"}]
    report = {
        "directory": str(lean),
        "files": len(files),
        "declarations": len(rows),
        "public": len(public),
        "private": len(private),
        "theorems": len(theorems),
        "names": [r["name"] for r in rows],
    }
    json.dump(report, sys.stdout, indent=2)
    sys.stdout.write("\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
