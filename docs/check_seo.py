#!/usr/bin/env python3
"""Validate SEO metadata on every generated HTML page."""

from __future__ import annotations

import sys
import tempfile
from html.parser import HTMLParser
from pathlib import Path
from urllib.parse import unquote, urlsplit

BASE_URL = "https://docs.dseams.info/"
REQUIRED_OG = {"og:title", "og:description", "og:image", "og:url", "og:type"}


class MetadataParser(HTMLParser):
    def __init__(self) -> None:
        super().__init__(convert_charrefs=True)
        self.in_title = False
        self.title = []
        self.descriptions = []
        self.canonicals = []
        self.open_graph = {}

    def handle_starttag(self, tag, attrs):
        values = {key.lower(): value for key, value in attrs if value is not None}
        if tag.lower() == "title":
            self.in_title = True
        elif tag.lower() == "meta":
            content = values.get("content", "").strip()
            if values.get("name", "").lower() == "description":
                self.descriptions.append(content)
            prop = values.get("property", "").lower()
            if prop.startswith("og:"):
                self.open_graph.setdefault(prop, []).append(content)
        elif tag.lower() == "link":
            rel = {item.lower() for item in values.get("rel", "").split()}
            if "canonical" in rel:
                self.canonicals.append(values.get("href", "").strip())

    def handle_endtag(self, tag):
        if tag.lower() == "title":
            self.in_title = False

    def handle_data(self, data):
        if self.in_title:
            self.title.append(data)


def _canonical_file(url: str, root: Path) -> Path | None:
    parsed = urlsplit(url)
    expected = urlsplit(BASE_URL)
    if (
        parsed.scheme != "https"
        or parsed.netloc != expected.netloc
        or parsed.query
        or parsed.fragment
    ):
        return None
    path = unquote(parsed.path)
    if path.endswith("/"):
        path += "index.html"
    relative = Path(path.lstrip("/"))
    if relative.is_absolute() or ".." in relative.parts:
        return None
    return root / relative


def page_problems(path: Path, root: Path) -> list[str]:
    parser = MetadataParser()
    parser.feed(path.read_text(encoding="utf-8", errors="replace"))
    label = path.relative_to(root).as_posix()
    problems = []
    if not "".join(parser.title).strip():
        problems.append(f"{label}: missing non-empty title")
    if not any(parser.descriptions):
        problems.append(f"{label}: missing non-empty meta description")
    if len(parser.canonicals) != 1:
        problems.append(f"{label}: expected one canonical URL, found {len(parser.canonicals)}")
    elif _canonical_file(parser.canonicals[0], root) != path:
        problems.append(f"{label}: canonical URL does not resolve to this page")
    for prop in sorted(REQUIRED_OG):
        values = parser.open_graph.get(prop, [])
        if len(values) != 1 or not values[0]:
            problems.append(f"{label}: expected one non-empty {prop}, found {len(values)}")
    if (
        len(parser.open_graph.get("og:url", [])) == 1
        and _canonical_file(parser.open_graph["og:url"][0], root) != path
    ):
        problems.append(f"{label}: og:url does not resolve to this page")
    return problems


def check_tree(root: Path) -> list[str]:
    pages = sorted(root.rglob("*.html")) if root.is_dir() else []
    if not pages:
        return [f"{root}: no built HTML pages found"]
    problems = []
    for page in pages:
        problems.extend(page_problems(page, root))
    return problems


def self_test() -> list[str]:
    with tempfile.TemporaryDirectory() as tmp:
        root = Path(tmp)
        good = root / "guide" / "index.html"
        good.parent.mkdir()
        good.write_text(
            """<title>Guide</title><meta name="description" content="Guide page">
<link rel="canonical" href="https://docs.dseams.info/guide/">
<meta property="og:title" content="Guide"><meta property="og:description" content="Guide page">
<meta property="og:image" content="https://dseams.info/images/og-card.png">
<meta property="og:url" content="https://docs.dseams.info/guide/">
<meta property="og:type" content="website">""",
            encoding="utf-8",
        )
        if page_problems(good, root):
            return ["self-test rejected valid metadata"]
        bad = root / "bad.html"
        bad.write_text("<html><head><title></title></head></html>", encoding="utf-8")
        found = page_problems(bad, root)
        expected = 3 + len(REQUIRED_OG)
        if len(found) != expected:
            return [f"self-test expected {expected} failures, found {len(found)}"]
    return []


def main(argv: list[str]) -> int:
    if argv[1:] == ["--self-test"]:
        problems = self_test()
    elif len(argv) == 2:
        problems = check_tree(Path(argv[1]))
    else:
        print("usage: docs/check_seo.py BUILD_DIR | --self-test", file=sys.stderr)
        return 2
    if problems:
        print("\n".join(problems), file=sys.stderr)
        return 1
    print("SEO metadata check ok")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
