#!/usr/bin/env python3
"""Reject maintainer process text and parse documented install claims offline."""

from __future__ import annotations

import html.parser
import re
import sys
from pathlib import Path

DOCS_DIR = Path(__file__).resolve().parent
if str(DOCS_DIR) not in sys.path:
    sys.path.insert(0, str(DOCS_DIR))

from registry_claims import ClaimParseError, claims

ROOT = Path(__file__).resolve().parents[1]

PROCESS_PATTERNS = (
    ("video claim file", re.compile(r"\bembed\.txt\b", re.I)),
    ("caption production note", re.compile(
        r"Captions are burned into the upload|Transcript:\s*docs/videos/|"
        r"docs/videos/[A-Za-z0-9_-]+/captions\.vtt",
        re.I,
    )),
    ("audio encoder flag", re.compile(r"\baac_coder\b", re.I)),
    ("repository policy", re.compile(r"\bstays? in the tree\b|\bout of the repository\b", re.I)),
    ("private cluster", re.compile(r"\bElja\b")),
    ("private cluster files", re.compile(r"\belja-(?:crt|sysroot)\b|elja_gpu_job\.sh", re.I)),
    ("job identifier", re.compile(r"\bjob [0-9]{5,}\b", re.I)),
    ("base versus tip", re.compile(r"\bBase [0-9a-f]{7,}\s+vs\b", re.I)),
    ("molecule rendering command", re.compile(r"\brender_molecules\.py\b|--out\s+docs/", re.I)),
    ("missing-player process", re.compile(r"\brenders? no player\b", re.I)),
    ("documentation test command", re.compile(r"\bmeson\s+test\b|\bdocs_[a-z0-9_]+\b", re.I)),
    ("continuous-integration process", re.compile(
        r"\bCI (?:executes|runs|mode)\b|\bin CI\b|\bcontinuous integration\b",
        re.I,
    )),
    ("workflow process", re.compile(r"\bSnakemake\b|\btowncrier\b|changelog\.d", re.I)),
    ("repository provenance", re.compile(
        r"\bpinned revision\b|\bsource_commit\b|"
        r"\b(?:commit|revision|recorded at|retimed at)\s+[`~]?[0-9a-f]{7,40}",
        re.I,
    )),
    ("internal source path", re.compile(r"docs/orgmode/contributing/repositories\.org", re.I)),
    ("internal run identifier", re.compile(
        r"\b(?:9fcaab36|cad887b2|e943930a|ac510645|18588dac|9896a478|"
        r"848426f9|caeb71f8|1ce3be3b|2ed49e87|6ec9c283|162cc834|"
        r"4828e33e|323bfd3872e726b3fa384a67c6b1bedf627ce6a0|"
        r"c7f71310f27c5ce62d11e97fd82f9757fa4e8628|"
        r"59f0546323e9e39ab7075b65969bbfd9ed4280bb|7f310225)\b|"
        r"\bis this sitting\b|\bpool-audit\b|\bHaoZeke/seams-core\b|"
        r"\bbench/h2h-[A-Za-z0-9_-]+\b",
        re.I,
    )),
    ("skip-as-success process", re.compile(
        r"\bexits? 77\b|--small\s+is\s+the\s+CI\s+mode|\breturns? success\b",
        re.I,
    )),
    ("documentation-production caveat", re.compile(
        r"\brunnable copy\b|\btest exits 77\b|\bprivate paper repository\b|"
        r"docs/records/README\.md|\bbeing uploaded\b|\bregistration is coming\b",
        re.I,
    )),
    ("binding publication process", re.compile(
        r"\bbeing published\b|\bregistry install\b.*\bcoming\b|"
        r"\bnot available there yet\b|"
        r"\bwith its API reference, tutorials,\s*tests and release notes\b|"
        r"\bthat this repository builds and installs\b|"
        r"\bscripts/install_pydseamslib\.sh\b",
        re.I,
    )),
    ("benchmark provenance", re.compile(
        r"\bbenchmarks/|\bdocs/(?:records|data)/|\bstored in git\b|"
        r"\braw lines\b|\bSource file\s+``|\bRendered from\b",
        re.I,
    )),
    ("registry publication status", re.compile(r"\bregistry publication\b", re.I)),
    ("unpublished LuaRocks package", re.compile(
        r"\bLuaRocks\b[^\n]{0,80}\bdseams\b",
        re.I,
    )),
    ("unpublished fpm package", re.compile(
        r"\bfpm\b(?:(?!\bGit\s+dependency\b)[^\n]){0,80}\bdseams\b",
        re.I,
    )),
    ("build status badge", re.compile(
        r"\bBuild Status\b|\bbuilt with nix\b|\bCheck our build status\b",
        re.I,
    )),
)

README_PATTERNS = tuple(
    item for item in PROCESS_PATTERNS if item[0] != "documentation test command"
)


def normalised_text(text: str) -> str:
    """Collapse horizontal whitespace before the reader-text patterns run.

    Sphinx splits a literal into one span per word and the visible-text
    join puts the span gaps back as several spaces. ``meson   test`` is
    the same command as ``meson test``. Newlines stay, so a pattern that
    is written to stop at the end of a line still does.
    """
    return re.sub(r"[^\S\n]+", " ", text)

SKIP_PARTS = {
    "science/figures",
    "videos/embeds",
}
SKIP_FILES = {
    "explanation/cutoff-rings-tables.rst",
    "explanation/performance-overview-tables.rst",
    "explanation/thread-scaling-tables.rst",
    "howto/pool-audit-v5-tables.rst",
    "reference/release-2.14-notes.rst",
    "videos/players.rst",
}


class VisibleHTML(html.parser.HTMLParser):
    def __init__(self):
        super().__init__()
        self.hidden = 0
        self.parts = []

    def handle_starttag(self, tag, attrs):
        if tag in {"script", "style", "template"}:
            self.hidden += 1

    def handle_endtag(self, tag):
        if tag in {"script", "style", "template"} and self.hidden:
            self.hidden -= 1

    def handle_data(self, data):
        if not self.hidden:
            self.parts.append(data)


def source_pages(root: Path | None = None):
    root = ROOT if root is None else root
    source = root / "docs" / "source"
    if source.is_dir():
        for path in sorted(source.rglob("*")):
            if path.suffix.lower() not in {".rst", ".md"}:
                continue
            rel = path.relative_to(source).as_posix()
            if rel in SKIP_FILES or any(rel == part or rel.startswith(part + "/") for part in SKIP_PARTS):
                continue
            yield path, rel
    llms = root / "llms.txt"
    if llms.is_file():
        yield llms, "llms.txt"
    readme = root / "README.md"
    if readme.is_file():
        yield readme, "README.md"


def html_pages(build):
    players = build / "videos" / "players.html"
    if players.exists():
        yield players, "videos/players.html", "players.html must not be published"
    for path in sorted(build.rglob("*.html")):
        rel = path.relative_to(build).as_posix()
        parser = VisibleHTML()
        parser.feed(path.read_text(encoding="utf-8", errors="replace"))
        yield path, rel, " ".join(parser.parts)


def page_findings(rel: str, text: str, scanned_text: str) -> tuple[list[str], int]:
    """Match one page. ``scanned_text`` is what the patterns see.

    ``main`` passes ``normalised_text(text)``. The claim parser still sees
    the original lines.
    """
    if text == "players.html must not be published":
        return [f"html/{rel}: internal fragment was published"], 0
    patterns = README_PATTERNS if rel == "README.md" else PROCESS_PATTERNS
    errors = []
    for label, pattern in patterns:
        for match in pattern.finditer(scanned_text):
            line = scanned_text.count("\n", 0, match.start()) + 1
            errors.append(f"{rel}:{line}: {label}: {match.group(0)}")
    registry_claims = 0
    for number, line in enumerate(text.splitlines(), 1):
        try:
            for _manager, _package, _version in claims(line):
                registry_claims += 1
        except ClaimParseError as error:
            errors.append(f"{rel}:{number}: malformed install claim: {error}")
    return errors, registry_claims


def main(root: Path | None = None) -> int:
    root = ROOT if root is None else root
    errors = []
    registry_claims = 0
    source_count = 0
    html_count = 0
    scanned = []
    for path, rel in source_pages(root):
        source_count += 1
        scanned.append((path, rel, path.read_text(encoding="utf-8", errors="replace")))
    build = root / "docs" / "build"
    if build.is_dir():
        for path, rel, text in html_pages(build):
            html_count += 1
            scanned.append((path, "html/" + rel, text))
    for path, rel, text in scanned:
        scanned_text = normalised_text(text)
        found, claims_here = page_findings(rel, text, scanned_text)
        errors.extend(found)
        registry_claims += claims_here
    if errors:
        print("\n".join(errors), file=sys.stderr)
        return 1
    print(
        f"reader text ok ({source_count} exported sources, {html_count} HTML pages; "
        f"{registry_claims} registry claims parsed)"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
