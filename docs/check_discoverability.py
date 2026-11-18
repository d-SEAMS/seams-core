#!/usr/bin/env python3
# Copyright (c) d-SEAMS developers
"""Check the discoverability files and, after a docs build, the HTML.

Source files (always):

- CITATION.cff and codemeta.json agree on version, license, and the paper
- context7.json names the manual folders and skips the vendored trees
- skills/d-seams/SKILL.md matches the Agent Skills frontmatter rules
- llms.txt matches the llms.txt section order

Built site (docs/build, the default run):

- sitemap.xml is sitemap XML and every URL has a lastmod
- robots.txt names the sitemap and allows the site
- 404.html is noindex and is not in the sitemap
- each sitemap HTML page has one title, description, canonical URL,
  Open Graph and Twitter card, JSON-LD, a Markdown twin, and stable
  heading anchors
- the social card is a 1200x630 PNG
- llms.txt and llms-full.txt were written by the build

Usage:
  docs/check_discoverability.py --self-test
  docs/check_discoverability.py
"""

import ast
import json
import re
import sys
import tempfile
from html.parser import HTMLParser
from pathlib import Path
from urllib.parse import urlsplit
from xml.etree import ElementTree as ET

ROOT = Path(__file__).resolve().parents[1]
DOCS = ROOT / "docs"
BUILD = DOCS / "build"
sys.path.insert(0, str(DOCS))

import discoverability as pub  # noqa: E402

SITE = pub.SITE
DOI = pub.DOI
HOST = "docs.dseams.info"
SITEMAP_NS = pub.SITEMAP_NS
LASTMOD = re.compile(
    r"^\d{4}-\d{2}-\d{2}"
    r"(?:T\d{2}:\d{2}:\d{2}(?:\.\d+)?(?:Z|[+-]\d{2}:\d{2}))?$"
)
LINK = re.compile(r"^- \[([^\]]+)\]\((https?://[^)\s]+)\)(?::\s+\S.*)?$")
HEADING = re.compile(r"^#{1,6} ")
H1 = re.compile(r"^# .+")
H2 = re.compile(r"^## .+")


def _version_problems():
    version = pub.project_version()
    problems = []
    if version != "2.14.0":
        problems.append("meson project version is {}, expected 2.14.0".format(version))
    conf = (DOCS / "source" / "conf.py").read_text(encoding="utf-8")
    if 'version = "{}"'.format(version) not in conf:
        problems.append("conf.py version does not match meson")
    return problems


def _citation_text():
    return (ROOT / "CITATION.cff").read_text(encoding="utf-8")


def citation_problems():
    text = _citation_text()
    version = pub.project_version()
    problems = []
    required = [
        "cff-version: 1.2.0",
        'version: "{}"'.format(version),
        "license: MIT",
        'repository-code: "https://github.com/d-SEAMS/seams-core"',
        'doi: "{}"'.format(DOI),
        "preferred-citation:",
        "journal: \"Journal of Chemical Information and Modeling\"",
        "year: 2020",
        "family-names: Goswami",
        "given-names: Rohit",
        "https://orcid.org/0000-0002-2393-8056",
        "date-released:",
    ]
    for needle in required:
        if needle not in text:
            problems.append("CITATION.cff lacks " + needle)
    released = re.search(r'^date-released:\s*"(\d{4}-\d{2}-\d{2})"', text, re.M)
    if not released:
        problems.append("CITATION.cff date-released is not YYYY-MM-DD")
    return problems


def codemeta_problems():
    path = ROOT / "codemeta.json"
    problems = []
    try:
        data = json.loads(path.read_text(encoding="utf-8"))
    except json.JSONDecodeError as exc:
        return ["codemeta.json is not JSON ({})".format(exc)]
    version = pub.project_version()
    if data.get("@type") != "SoftwareSourceCode":
        problems.append("codemeta @type is not SoftwareSourceCode")
    if data.get("version") != version:
        problems.append("codemeta version does not match meson")
    if data.get("license") != pub.LICENSE:
        problems.append("codemeta license is not the MIT SPDX URL")
    if data.get("codeRepository") != pub.REPO:
        problems.append("codemeta codeRepository is not the seams-core URL")
    if DOI not in json.dumps(data):
        problems.append("codemeta does not cite {}".format(DOI))
    languages = data.get("programmingLanguage")
    if not isinstance(languages, list) or "C++" not in languages or "C" not in languages:
        problems.append("codemeta programmingLanguage must list C++ and C")
    authors = data.get("author")
    names = []
    if isinstance(authors, list):
        for person in authors:
            if isinstance(person, dict):
                names.append(
                    "{} {}".format(person.get("givenName", ""), person.get("familyName", ""))
                )
    if "Rohit Goswami" not in names:
        problems.append("codemeta authors omit Rohit Goswami")
    released = re.search(
        r'^date-released:\s*"(\d{4}-\d{2}-\d{2})"', _citation_text(), re.M
    )
    if released and data.get("datePublished") != released.group(1):
        problems.append("codemeta datePublished does not match CITATION.cff")
    return problems


def context7_problems():
    path = ROOT / "context7.json"
    try:
        data = json.loads(path.read_text(encoding="utf-8"))
    except json.JSONDecodeError as exc:
        return ["context7.json is not JSON ({})".format(exc)]
    problems = []
    if data.get("$schema") != "https://context7.com/schema/context7.json":
        problems.append("context7.json schema is not the Context7 schema")
    title = data.get("projectTitle", "")
    if not isinstance(title, str) or not 1 <= len(title) <= 100:
        problems.append("context7 projectTitle must be 1 to 100 characters")
    description = data.get("description", "")
    if not isinstance(description, str) or not 1 <= len(description) <= 500:
        problems.append("context7 description must be 1 to 500 characters")
    elif description.lower().startswith(("ignore ", "you are ", "as an ai")):
        problems.append("context7 description is an instruction, not a description")
    folders = data.get("folders")
    if not isinstance(folders, list) or not folders:
        problems.append("context7 folders is empty")
        folders = []
    for folder in folders:
        if not isinstance(folder, str) or not (ROOT / folder).is_dir():
            problems.append("context7 folder is missing: {}".format(folder))
    for needed in ("docs/orgmode", "docs/api", "docs/tutorials", "skills"):
        if needed not in folders:
            problems.append("context7 folders omit " + needed)
    excluded = data.get("excludeFolders") or []
    for needed in ("subprojects", "third_party", "tests", "docs/source", "docs/build"):
        if needed not in excluded:
            problems.append("context7 excludeFolders omit " + needed)
    rules = data.get("rules") or []
    if not isinstance(rules, list) or not rules:
        problems.append("context7 rules is empty")
    for rule in rules:
        if not isinstance(rule, str) or not 5 <= len(rule) <= 200:
            problems.append("context7 rule length is outside 5 to 200: {}".format(rule))
    return problems


def skill_problems():
    path = ROOT / "skills" / "d-seams" / "SKILL.md"
    problems = []
    if not path.is_file():
        return ["skills/d-seams/SKILL.md is missing"]
    text = path.read_text(encoding="utf-8")
    if not text.startswith("---\n"):
        problems.append("SKILL.md does not open with frontmatter")
        return problems
    end = text.find("\n---\n", 3)
    if end < 0:
        problems.append("SKILL.md frontmatter does not close")
        return problems
    front = text[4:end]
    body = text[end + 5 :]
    name = re.search(r"^name:\s*([A-Za-z0-9-]+)\s*$", front, re.M)
    if name is None or name.group(1) != "d-seams" or path.parent.name != "d-seams":
        problems.append("SKILL.md name must be d-seams and match the directory")
    elif "--" in name.group(1) or name.group(1).startswith("-") or name.group(1).endswith("-"):
        problems.append("SKILL.md name has a bad hyphen")
    desc = re.search(r'^description:\s*"([^"]*)"', front, re.M)
    if desc is None or not 1 <= len(desc.group(1)) <= 1024:
        problems.append("SKILL.md description must be a quoted string of 1 to 1024 characters")
    elif "seams" not in desc.group(1).lower() or "when" not in desc.group(1).lower():
        problems.append("SKILL.md description should say what it covers and when to use it")
    if 'version: "{}"'.format(pub.project_version()) not in front:
        problems.append("SKILL.md version does not match meson")
    if "license: MIT" not in front:
        problems.append("SKILL.md license is not MIT")
    for needle in ("seams chill-plus", "pydseamslib", "seams.h", "3.5"):
        if needle not in body and needle not in text:
            problems.append("SKILL.md does not mention " + needle)
    if len(body.strip()) < 200:
        problems.append("SKILL.md body is too short to brief an agent")
    return problems


def agents_problems():
    path = ROOT / "AGENTS.md"
    if not path.is_file():
        return ["AGENTS.md is missing"]
    text = path.read_text(encoding="utf-8")
    problems = []
    for needle in (
        "skills/d-seams/SKILL.md",
        "llms.txt",
        "context7.json",
        "Do not weaken, skip, or",
    ):
        if needle not in text:
            problems.append("AGENTS.md lacks " + needle)
    return problems


def readme_problems():
    text = (ROOT / "README.md").read_text(encoding="utf-8")
    problems = []
    for needle in (
        "https://anaconda.org/conda-forge/libseams",
        "https://pypi.org/project/pydseamslib/",
        "https://docs.dseams.info",
        "https://doi.org/" + DOI,
    ):
        if needle not in text:
            problems.append("README lacks " + needle)
    if not text.startswith("# d-SEAMS\n"):
        problems.append("README does not open with the project name")
    return problems


def llms_problems(text, label="llms.txt"):
    """The llms.txt section order. Returns a list of problems."""
    if text.startswith("\ufeff"):
        text = text[1:]
    lines = text.splitlines()
    if not lines or not H1.match(lines[0]) or lines[0].startswith("##"):
        return ["{} must open with one H1".format(label)]
    problems = []
    phase = "quote"
    saw_quote = False
    saw_h2 = False
    for number, line in enumerate(lines[1:], start=2):
        if phase == "quote":
            if line.startswith(">"):
                saw_quote = True
                continue
            if line == "":
                if saw_quote:
                    phase = "free"
                continue
            problems.append("{}:{} blockquote is missing after the H1".format(label, number))
            phase = "free"
        if HEADING.match(line) and not H2.match(line):
            problems.append("{}:{} heading is not an H2 section".format(label, number))
            continue
        if H2.match(line):
            phase = "files"
            saw_h2 = True
            continue
        if phase == "free":
            continue
        if line == "":
            continue
        if not LINK.match(line):
            problems.append("{}:{} file list entry is not a markdown link".format(label, number))
    if not saw_quote:
        problems.append("{} has no summary blockquote".format(label))
    if not saw_h2:
        problems.append("{} has no H2 file list".format(label))
    return problems


def _local(tag):
    if tag.startswith("{"):
        return tag.partition("}")[2]
    return tag


def sitemap_records(path):
    """Return (problems, [(loc, lastmod), ...])."""
    try:
        root = ET.parse(path).getroot()
    except ET.ParseError as exc:
        return ["sitemap.xml is not XML ({})".format(exc)], []
    if _local(root.tag) != "urlset":
        return ["sitemap.xml root is not urlset"], []
    problems = []
    rows = []
    for child in list(root):
        if _local(child.tag) != "url":
            continue
        loc = ""
        last = ""
        for field in list(child):
            if _local(field.tag) == "loc" and field.text:
                loc = field.text.strip()
            elif _local(field.tag) == "lastmod" and field.text:
                last = field.text.strip()
        if not loc:
            problems.append("sitemap has a url with no loc")
            continue
        parsed = urlsplit(loc)
        if parsed.scheme != "https" or parsed.hostname != HOST:
            problems.append("sitemap loc is off-site: " + loc)
        if not last or not LASTMOD.match(last):
            problems.append("sitemap loc has no lastmod: " + loc)
        rows.append((loc, last))
    if not rows:
        problems.append("sitemap has no URLs")
    return problems, rows


class _Page(HTMLParser):
    def __init__(self):
        super().__init__(convert_charrefs=True)
        self.title = []
        self.in_title = False
        self.metas = []
        self.links = []
        self.scripts = []
        self.in_json = False
        self.json_buf = []
        self.lang = ""
        self.article = 0
        self.sections = []
        self.headings = []
        self._heading = None
        self._heading_buf = []

    def handle_starttag(self, tag, attrs):
        attr = {key.lower(): value for key, value in attrs}
        if tag == "html":
            self.lang = attr.get("lang", "")
        elif tag == "title":
            self.in_title = True
        elif tag == "meta":
            self.metas.append(attr)
        elif tag == "link":
            self.links.append(attr)
        elif tag == "script" and attr.get("type") == "application/ld+json":
            self.in_json = True
            self.json_buf = []
        elif tag == "article":
            self.article += 1
        elif tag == "section":
            self.sections.append(attr.get("id", ""))
        elif self.article and tag in ("h1", "h2", "h3", "h4", "h5", "h6") and self._heading is None:
            ident = attr.get("id") or (self.sections[-1] if self.sections else "")
            self._heading = (tag, ident)
            self._heading_buf = []

    def handle_endtag(self, tag):
        if tag == "title":
            self.in_title = False
        elif tag == "script" and self.in_json:
            self.scripts.append("".join(self.json_buf))
            self.in_json = False
        elif tag == "article" and self.article:
            self.article -= 1
        elif tag == "section" and self.sections:
            self.sections.pop()
        elif self._heading and tag == self._heading[0]:
            self.headings.append((self._heading[0], self._heading[1], "".join(self._heading_buf).strip()))
            self._heading = None

    def handle_data(self, data):
        if self.in_title:
            self.title.append(data)
        if self.in_json:
            self.json_buf.append(data)
        if self._heading:
            self._heading_buf.append(data)


def _meta(page, name):
    found = []
    for attr in page.metas:
        if attr.get("name", "").lower() == name and "content" in attr:
            found.append(attr["content"])
    return found


def _prop(page, name):
    found = []
    for attr in page.metas:
        if attr.get("property", "").lower() == name and "content" in attr:
            found.append(attr["content"])
    return found


def _links(page, rel, type_=None):
    found = []
    for attr in page.links:
        if attr.get("rel", "").lower() != rel:
            continue
        if type_ is not None and attr.get("type", "").lower() != type_:
            continue
        if attr.get("href"):
            found.append(attr["href"])
    return found


def _types(node, found):
    if isinstance(node, dict):
        kind = node.get("@type")
        if isinstance(kind, str):
            found.setdefault(kind, []).append(node)
        elif isinstance(kind, list):
            for item in kind:
                if isinstance(item, str):
                    found.setdefault(item, []).append(node)
        for value in node.values():
            _types(value, found)
    elif isinstance(node, list):
        for item in node:
            _types(item, found)


def _software_problems(node, label):
    problems = []
    version = node.get("version") or node.get("softwareVersion")
    if version != pub.project_version():
        problems.append("{} version is not the project version".format(label))
    if node.get("license") != pub.LICENSE:
        problems.append("{} license is missing".format(label))
    languages = node.get("programmingLanguage")
    if not isinstance(languages, list) or len(languages) < 2:
        problems.append("{} programmingLanguage is not a list".format(label))
    citation = node.get("citation")
    blob = json.dumps(citation) if citation is not None else ""
    if DOI not in blob:
        problems.append("{} citation lacks the paper DOI".format(label))
    return problems


# Doxygen writes api-cpp/ after Sphinx. Those pages have no Markdown twin or
# manual JSON-LD. check_seo.py covers their metadata.
DOXYGEN_DIR = "api-cpp/"


def _route_files():
    """Map each public route (no ``.html``, lower case) to its built page."""
    routes = {}
    for path in BUILD.rglob("*.html"):
        route = "/" + path.relative_to(BUILD).as_posix()
        if route.endswith("/index.html"):
            route = route[:-10] or "/"
        else:
            route = route[:-5]
        routes[route.lower()] = path
    return routes


def page_problems(html_path, md_path):
    parser = _Page()
    try:
        parser.feed(html_path.read_text(encoding="utf-8", errors="replace"))
    except Exception as exc:
        return ["{} did not parse ({})".format(html_path.name, exc)]
    problems = []
    title = " ".join("".join(parser.title).split())
    rel = html_path.relative_to(BUILD).as_posix()
    if parser.lang != "en":
        problems.append("{} lang is not en".format(rel))
    if not title or "d-SEAMS" not in title:
        problems.append("{} title is missing the project name".format(rel))
    twitter_title = _meta(parser, "twitter:title")
    if twitter_title != [title]:
        problems.append("{} twitter:title does not match the title".format(rel))
    descriptions = _meta(parser, "description")
    if len(descriptions) != 1:
        problems.append("{} needs one meta description".format(rel))
    elif not 40 <= len(descriptions[0]) <= 180:
        problems.append("{} description length is {}".format(rel, len(descriptions[0])))
    canonical = _links(parser, "canonical")
    expected = SITE + "/" + rel
    if canonical != [expected]:
        problems.append("{} canonical is {}".format(rel, canonical))
    for prop in ("og:title", "og:description", "og:url", "og:type", "og:image"):
        if len(_prop(parser, prop)) != 1:
            problems.append("{} needs one {}".format(rel, prop))
    if _prop(parser, "og:url") not in ([], [expected]):
        # og:url is the canonical page URL.
        if _prop(parser, "og:url") != [expected]:
            problems.append("{} og:url is not the canonical URL".format(rel))
    if _prop(parser, "og:image") != [pub.OG_URL]:
        problems.append("{} og:image is not the social card".format(rel))
    if _meta(parser, "twitter:card") != ["summary_large_image"]:
        problems.append("{} twitter card is not summary_large_image".format(rel))
    if _meta(parser, "twitter:image") != [pub.OG_URL]:
        problems.append("{} twitter:image is not the social card".format(rel))
    if _meta(parser, "twitter:description") != descriptions:
        problems.append("{} twitter description does not match".format(rel))
    robots = " ".join(_meta(parser, "robots")).lower()
    if "noindex" in robots:
        problems.append("{} is noindex".format(rel))
    markdown = _links(parser, "alternate", "text/markdown")
    md_url = SITE + "/" + rel[:-5] + ".md" if rel.endswith(".html") else ""
    if markdown != [md_url]:
        problems.append("{} markdown alternate is {}".format(rel, markdown))
    if _links(parser, "describedby") != [SITE + "/llms.txt"]:
        problems.append("{} does not point at llms.txt".format(rel))
    if len(parser.scripts) != 1:
        problems.append("{} needs one JSON-LD script".format(rel))
    else:
        try:
            data = json.loads(parser.scripts[0])
        except json.JSONDecodeError as exc:
            problems.append("{} JSON-LD did not parse ({})".format(rel, exc))
            data = None
        if data is not None:
            found = {}
            _types(data, found)
            for kind in ("SoftwareSourceCode", "SoftwareApplication", "WebPage"):
                if kind not in found:
                    problems.append("{} JSON-LD lacks {}".format(rel, kind))
            for kind in ("SoftwareSourceCode", "SoftwareApplication"):
                for node in found.get(kind, []):
                    problems.extend(
                        "{} {}".format(rel, item) for item in _software_problems(node, kind)
                    )
    if not md_path.is_file():
        problems.append("{} has no Markdown twin".format(rel))
        return problems
    markdown_text = md_path.read_text(encoding="utf-8", errors="replace")
    seen = set()
    for _level, ident, _text in parser.headings:
        if not ident or not pub.ID_OK.match(ident):
            problems.append("{} heading has a non-slug id {!r}".format(rel, ident))
            continue
        if ident in seen:
            problems.append("{} repeats heading id {}".format(rel, ident))
        seen.add(ident)
        if "{#" + ident + "}" not in markdown_text:
            problems.append("{} markdown lacks anchor {}".format(rel, ident))
    if not parser.headings:
        problems.append("{} article has no headings".format(rel))
    return problems


def build_problems():
    problems = []
    if not BUILD.is_dir():
        return ["docs/build is missing. Run pixi run -e docs docbld"]
    robots = BUILD / "robots.txt"
    if not robots.is_file():
        problems.append("robots.txt is missing")
    else:
        text = robots.read_text(encoding="utf-8")
        if "Sitemap: {}/sitemap.xml".format(SITE) not in text:
            problems.append("robots.txt does not name the sitemap")
        if "Allow: /" not in text or re.search(r"(?m)^Disallow:\s*/\s*$", text):
            problems.append("robots.txt does not allow the site")
        if "search=yes" not in text or "ai-input=yes" not in text:
            problems.append("robots.txt does not allow search and AI input")
    missing = BUILD / "404.html"
    if not missing.is_file():
        problems.append("404.html is missing")
    else:
        page = missing.read_text(encoding="utf-8", errors="replace")
        if "noindex" not in page.lower():
            problems.append("404.html is not noindex")
        if "Page not found" not in page or "{}/index.html".format(SITE) not in page:
            problems.append("404.html does not link back to the manual")
    image = BUILD / pub.OG_PATH
    if not image.is_file():
        problems.append("social card is missing")
    else:
        size = pub.png_size(image.read_bytes())
        if size != (1200, 630):
            problems.append("social card is {} not 1200x630".format(size))
    sitemap = BUILD / "sitemap.xml"
    if not sitemap.is_file():
        return problems + ["sitemap.xml is missing"]
    site_problems, rows = sitemap_records(sitemap)
    problems.extend(site_problems)
    locs = [loc for loc, _last in rows]
    if any(loc.rstrip("/").endswith("/404.html") for loc in locs):
        problems.append("sitemap lists 404.html")
    for banned in ("/search.html", "/genindex.html"):
        if any(loc.endswith(banned) for loc in locs):
            problems.append("sitemap lists " + banned)
    titles = {}
    descriptions = {}
    files = _route_files()
    for loc, _last in rows:
        html_path = files.get(urlsplit(loc).path.lower())
        if html_path is None:
            problems.append("sitemap loc has no file: " + loc)
            continue
        rel = html_path.relative_to(BUILD).as_posix()
        if rel.startswith(DOXYGEN_DIR):
            continue
        md_path = BUILD / (rel[:-5] + ".md")
        got = page_problems(html_path, md_path)
        problems.extend(got)
        if html_path.is_file():
            parser = _Page()
            parser.feed(html_path.read_text(encoding="utf-8", errors="replace"))
            title = " ".join("".join(parser.title).split())
            titles.setdefault(title, []).append(rel)
            desc = _meta(parser, "description")
            if desc:
                descriptions.setdefault(desc[0], []).append(rel)
    for title, pages in titles.items():
        if title and len(pages) > 1:
            problems.append("title {!r} is used by {}".format(title, ", ".join(pages)))
    for desc, pages in descriptions.items():
        if len(pages) > 8:
            problems.append("one description is reused {} times".format(len(pages)))
    llms_path = BUILD / "llms.txt"
    if not llms_path.is_file():
        problems.append("built llms.txt is missing")
    else:
        built = llms_path.read_text(encoding="utf-8")
        problems.extend(llms_problems(built, "built llms.txt"))
        fresh = pub._llms().render()
        if built != fresh:
            problems.append("built llms.txt is not the generator output")
        for match in re.finditer(r"\((https://docs\.dseams\.info/[^)\s]+)\)", built):
            url = match.group(1)
            target = urlsplit(url).path.lstrip("/")
            if target.endswith(".md") and not (BUILD / target).is_file():
                problems.append("llms.txt links a missing file: " + url)
    full = BUILD / "llms-full.txt"
    if not full.is_file():
        problems.append("llms-full.txt is missing")
    else:
        problems.extend(pub.full_manual_problems(
            full.read_text(encoding="utf-8", errors="replace")))
    return problems


def fixture_problems():
    """Parsers on synthetic inputs, so the check is tested without a book build."""
    problems = []
    good = (
        "# Example\n\n"
        "> A short summary.\n\n"
        "One paragraph of detail.\n\n"
        "## Docs\n\n"
        "- [Start](https://example.com/start.md): The first page\n\n"
        "## Optional\n\n"
        "- [Full](https://example.com/llms-full.txt): Everything\n"
    )
    if llms_problems(good, "fixture"):
        problems.append("valid llms.txt was rejected: " + "; ".join(llms_problems(good, "fixture")))
    bad = "# Example\n\n## Early\n\nNo quote.\n"
    if not llms_problems(bad, "fixture"):
        problems.append("llms.txt without a summary was accepted")
    stuffed = "# Example\n\n> Fine.\n\n## Docs\n\njust a sentence\n"
    if not any("file list" in item for item in llms_problems(stuffed, "fixture")):
        problems.append("a file list without a link was accepted")
    size = pub.png_size(pub.social_png())
    if size != (1200, 630):
        problems.append("social card generator returned {}".format(size))
    try:
        data = json.loads(pub.json_ld("quickstart", "Quickstart - d-SEAMS", "Install the engine."))
    except json.JSONDecodeError as exc:
        problems.append("json_ld did not parse ({})".format(exc))
        data = None
    if data is not None:
        found = {}
        _types(data, found)
        for kind in ("SoftwareSourceCode", "SoftwareApplication"):
            if kind not in found:
                problems.append("json_ld lacks " + kind)
            for node in found.get(kind, []):
                problems.extend(_software_problems(node, kind))
    with tempfile.TemporaryDirectory() as tmp:
        sitemap = Path(tmp) / "sitemap.xml"
        sitemap.write_text(
            "<?xml version='1.0' encoding='utf-8'?>\n"
            '<urlset xmlns="http://www.sitemaps.org/schemas/sitemap/0.9">\n'
            "<url><loc>https://docs.dseams.info/index.html</loc></url>\n"
            "<url><loc>https://docs.dseams.info/quickstart.html</loc>"
            "<lastmod>2026-11-15T00:00:00Z</lastmod></url>\n"
            "</urlset>\n",
            encoding="utf-8",
        )
        pub.stamp_sitemap(sitemap, lambda _loc: "2026-11-15T00:00:00+00:00")
        got, rows = sitemap_records(sitemap)
        if got or len(rows) != 2 or any(not last for _loc, last in rows):
            problems.append("stamp_sitemap did not fill lastmod ({})".format("; ".join(got)))
    clipped = pub.clip("word " * 80, 40)
    if clipped.endswith(" ") or len(clipped) > 40:
        problems.append("clip exceeded its limit")
    return problems


def source_problems():
    problems = []
    problems.extend(_version_problems())
    problems.extend(citation_problems())
    problems.extend(codemeta_problems())
    problems.extend(context7_problems())
    problems.extend(skill_problems())
    problems.extend(agents_problems())
    problems.extend(readme_problems())
    llms = ROOT / "llms.txt"
    if not llms.is_file():
        problems.append("llms.txt is missing")
    else:
        problems.extend(llms_problems(llms.read_text(encoding="utf-8"), "llms.txt"))
    return problems


def linkcheck_manual_problems():
    """The linkcheck stand-in is a real manual, and other URLs stay remote."""
    problems = []

    class Html:
        class builder:
            name = "html"

    if pub._on_linkcheck_uri(None, pub.SITE + "/llms.txt") is not None:
        problems.append("linkcheck rewrote llms.txt")
    if pub._on_linkcheck_uri(Html(), pub.SITE + "/llms-full.txt") is not None:
        problems.append("the HTML builder rewrote llms-full.txt")
    resolved = pub._on_linkcheck_uri(None, pub.SITE + "/llms-full.txt")
    if not resolved or not Path(resolved).is_file():
        return problems + ["linkcheck did not resolve llms-full.txt to a file"]
    text = Path(resolved).read_text(encoding="utf-8")
    problems.extend(pub.full_manual_problems(text))
    if "#+TITLE:" not in text:
        problems.append("linkcheck manual is not assembled from the Org sources")
    return problems


def self_test():
    problems = []
    problems.extend(source_problems())
    problems.extend(fixture_problems())
    problems.extend(linkcheck_manual_problems())
    # conf.py still parses. The book imports it; a syntax error fails the build.
    ast.parse((DOCS / "source" / "conf.py").read_text(encoding="utf-8"))
    return problems


def main(argv):
    if argv[1:] == ["--self-test"]:
        problems = self_test()
    elif not argv[1:]:
        problems = source_problems()
        problems.extend(build_problems())
    else:
        print("usage: docs/check_discoverability.py [--self-test]", file=sys.stderr)
        return 2
    if problems:
        print("\n".join(problems), file=sys.stderr)
        return 1
    print("discoverability check ok")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
