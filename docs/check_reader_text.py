#!/usr/bin/env python3
"""Reject maintainer process text and unrecorded package promises."""

from __future__ import annotations

import csv
import re
import shlex
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "docs" / "source"
VERSIONS = ROOT / "docs" / "data" / "published-versions.tsv"

PROCESS_PATTERNS = (
    ("video claim file", re.compile(r"\bembed\.txt\b", re.I)),
    ("audio encoder flag", re.compile(r"\baac_coder\b", re.I)),
    ("repository policy", re.compile(r"\bstays? in the tree\b|\bout of the repository\b", re.I)),
    ("private cluster", re.compile(r"\bElja\b")),
    ("job identifier", re.compile(r"\bjob [0-9]{5,}\b", re.I)),
    ("base versus tip", re.compile(r"\bBase [0-9a-f]{7,}\s+vs\b", re.I)),
    ("molecule rendering command", re.compile(r"\brender_molecules\.py\b|--out\s+docs/", re.I)),
    ("missing-player process", re.compile(r"\brenders? no player\b", re.I)),
    ("documentation test command", re.compile(r"\bmeson test\b|\bdocs_[a-z0-9_]+\b", re.I)),
    ("continuous-integration process", re.compile(r"\bCI\b|continuous integration", re.I)),
    ("workflow process", re.compile(r"\bSnakemake\b|\btowncrier\b|changelog\.d", re.I)),
    ("repository provenance", re.compile(
        r"\bpinned revision\b|\bsource_commit\b|"
        r"\b(?:commit|revision|recorded at|retimed at)\s+[`~]?[0-9a-f]{7,40}",
        re.I,
    )),
    ("documentation-production caveat", re.compile(
        r"\brunnable copy\b|\btest exits 77\b|\bprivate paper repository\b|"
        r"docs/records/README\.md|\bbeing uploaded\b|\bregistration is coming\b",
        re.I,
    )),
)

SKIP_PARTS = {
    "contributing",
    "science/figures",
    "videos/embeds",
}
SKIP_FILES = {
    "reference/release-2.14-notes.rst",
    "videos/players.rst",
}


def published_versions() -> set[tuple[str, str, str]]:
    with VERSIONS.open(newline="", encoding="utf-8") as handle:
        rows = csv.DictReader(handle, delimiter="\t")
        required = {"manager", "package", "version"}
        if set(rows.fieldnames or ()) != required:
            raise ValueError(f"{VERSIONS}: expected columns {sorted(required)}")
        return {
            (row["manager"].lower(), row["package"].lower(), row["version"])
            for row in rows
        }


def claims(line: str):
    plain = line.replace("~", "").replace("`", "")

    for match in re.finditer(r"\bpip(?:3)?\s+install(?:\s+-\S+)*\s+['\"]?([A-Za-z0-9_.-]+)(?:==([A-Za-z0-9_.+-]+))?", plain):
        yield "pip", match.group(1), match.group(2) or "*"

    match = re.search(r"\b(?:conda|mamba)\s+install\b(.*)", plain)
    if match:
        try:
            tokens = shlex.split(match.group(1))
        except ValueError:
            tokens = match.group(1).split()
        skip = False
        for token in tokens:
            if skip:
                skip = False
                continue
            if token in {"-c", "--channel"}:
                skip = True
                continue
            if token.startswith("-"):
                continue
            package, sep, version = token.partition("=")
            yield "conda", package, version if sep else "*"

    for match in re.finditer(r"\bluarocks\s+install\s+([A-Za-z0-9_.-]+)(?:\s+([A-Za-z0-9_.+-]+))?", plain):
        yield "luarocks", match.group(1), match.group(2) or "*"

    for call in re.finditer(r"\bPkg\.add\(([^)]*)\)", plain):
        body = call.group(1)
        if "url=" in body:
            continue
        name = re.search(r'name\s*=\s*["\']([^"\']+)', body)
        version = re.search(r'version\s*=\s*["\']([^"\']+)', body)
        if name:
            yield "julia", name.group(1), version.group(1) if version else "*"
        else:
            short = re.fullmatch(r'\s*["\']([^"\']+)["\']\s*', body)
            if short:
                yield "julia", short.group(1), "*"

    for match in re.finditer(r"\binstall\.packages\(\s*[\"']([^\"']+)[\"']", plain):
        yield "r", match.group(1), "*"


def pages():
    for path in sorted(SOURCE.rglob("*")):
        if path.suffix.lower() not in {".rst", ".md"}:
            continue
        rel = path.relative_to(SOURCE).as_posix()
        if rel in SKIP_FILES or any(rel == part or rel.startswith(part + "/") for part in SKIP_PARTS):
            continue
        yield path, rel
    llms = ROOT / "llms.txt"
    if llms.is_file():
        yield llms, "llms.txt"


def main() -> int:
    allowed = published_versions()
    errors = []
    count = 0
    for path, rel in pages():
        count += 1
        text = path.read_text(encoding="utf-8", errors="replace")
        for label, pattern in PROCESS_PATTERNS:
            for match in pattern.finditer(text):
                line = text.count("\n", 0, match.start()) + 1
                errors.append(f"{rel}:{line}: {label}: {match.group(0)}")
        for number, line in enumerate(text.splitlines(), 1):
            for manager, package, version in claims(line):
                key = (manager, package.lower(), version)
                if key not in allowed:
                    errors.append(
                        f"{rel}:{number}: unrecorded install claim: "
                        f"{manager} {package} {version}"
                    )
    if errors:
        print("\n".join(errors), file=sys.stderr)
        return 1
    print(f"reader text ok ({count} exported pages)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
