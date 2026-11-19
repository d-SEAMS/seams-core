#!/usr/bin/env python3
# Copyright (c) d-SEAMS developers
"""Publish the discoverability files for the d-SEAMS manual.

The HTML builder writes, for every page:

- a description, a canonical URL, Open Graph and Twitter tags
- schema.org JSON-LD for the software and the page
- a Markdown twin (``page.md`` beside ``page.html``)
- ``rel=alternate`` and ``rel=describedby`` links

When the HTML build finishes it also writes ``robots.txt``, ``404.html``,
``llms.txt``, ``llms-full.txt``, the social card, and a ``lastmod`` on
every sitemap URL. Sphinx already writes the sitemap. This module fills
the dates from git, because the generated ``.rst`` is not a tracked file
and ``sphinx-last-updated-by-git`` is not part of the docs environment.
"""

from __future__ import annotations

import importlib.util
import json
import re
import struct
import subprocess
import zlib
from pathlib import Path
from urllib.parse import urlsplit
from xml.etree import ElementTree as ET

ROOT = Path(__file__).resolve().parents[1]
SITE = "https://docs.dseams.info"
DOI = "10.1021/acs.jcim.0c00031"
DOI_URL = "https://doi.org/" + DOI
LICENSE = "https://spdx.org/licenses/MIT"
REPO = "https://github.com/d-SEAMS/seams-core"
OG_PATH = "_static/og-dseams.png"
OG_URL = SITE + "/" + OG_PATH
OG_ALT = "d-SEAMS. Deferred structural elucidation analysis for molecular simulations."
SITEMAP_NS = "http://www.sitemaps.org/schemas/sitemap/0.9"
ID_OK = re.compile(r"^[A-Za-z][A-Za-z0-9_-]*$")
_LLMS = None
_GIT_TIME = {}
_HEAD_TIME = None

# 5x7 glyphs, high bit on the left. The social card is drawn from these.
_FONT = {
    " ": (0, 0, 0, 0, 0, 0, 0),
    "A": (0b01110, 0b10001, 0b10001, 0b11111, 0b10001, 0b10001, 0b10001),
    "B": (0b11110, 0b10001, 0b10001, 0b11110, 0b10001, 0b10001, 0b11110),
    "C": (0b01111, 0b10000, 0b10000, 0b10000, 0b10000, 0b10000, 0b01111),
    "D": (0b11110, 0b10001, 0b10001, 0b10001, 0b10001, 0b10001, 0b11110),
    "E": (0b11111, 0b10000, 0b10000, 0b11110, 0b10000, 0b10000, 0b11111),
    "F": (0b11111, 0b10000, 0b10000, 0b11110, 0b10000, 0b10000, 0b10000),
    "G": (0b01111, 0b10000, 0b10000, 0b10011, 0b10001, 0b10001, 0b01111),
    "H": (0b10001, 0b10001, 0b10001, 0b11111, 0b10001, 0b10001, 0b10001),
    "I": (0b01110, 0b00100, 0b00100, 0b00100, 0b00100, 0b00100, 0b01110),
    "J": (0b00111, 0b00010, 0b00010, 0b00010, 0b10010, 0b10010, 0b01100),
    "K": (0b10001, 0b10010, 0b10100, 0b11000, 0b10100, 0b10010, 0b10001),
    "L": (0b10000, 0b10000, 0b10000, 0b10000, 0b10000, 0b10000, 0b11111),
    "M": (0b10001, 0b11011, 0b10101, 0b10101, 0b10001, 0b10001, 0b10001),
    "N": (0b10001, 0b11001, 0b10101, 0b10011, 0b10001, 0b10001, 0b10001),
    "O": (0b01110, 0b10001, 0b10001, 0b10001, 0b10001, 0b10001, 0b01110),
    "P": (0b11110, 0b10001, 0b10001, 0b11110, 0b10000, 0b10000, 0b10000),
    "Q": (0b01110, 0b10001, 0b10001, 0b10001, 0b10101, 0b10010, 0b01101),
    "R": (0b11110, 0b10001, 0b10001, 0b11110, 0b10100, 0b10010, 0b10001),
    "S": (0b01111, 0b10000, 0b10000, 0b01110, 0b00001, 0b00001, 0b11110),
    "T": (0b11111, 0b00100, 0b00100, 0b00100, 0b00100, 0b00100, 0b00100),
    "U": (0b10001, 0b10001, 0b10001, 0b10001, 0b10001, 0b10001, 0b01110),
    "V": (0b10001, 0b10001, 0b10001, 0b10001, 0b10001, 0b01010, 0b00100),
    "W": (0b10001, 0b10001, 0b10001, 0b10101, 0b10101, 0b10101, 0b01010),
    "X": (0b10001, 0b10001, 0b01010, 0b00100, 0b01010, 0b10001, 0b10001),
    "Y": (0b10001, 0b10001, 0b01010, 0b00100, 0b00100, 0b00100, 0b00100),
    "Z": (0b11111, 0b00001, 0b00010, 0b00100, 0b01000, 0b10000, 0b11111),
    "+": (0b00100, 0b00100, 0b11111, 0b00100, 0b00100, 0, 0),
    "-": (0, 0, 0, 0b11111, 0, 0, 0),
    ".": (0, 0, 0, 0, 0, 0b01100, 0b01100),
}


def project_version():
    text = (ROOT / "meson.build").read_text(encoding="utf-8")
    match = re.search(r"version:\s*'([0-9]+\.[0-9]+\.[0-9]+)'", text)
    if not match:
        raise RuntimeError("meson.build has no project version")
    return match.group(1)


def clip(text, limit=155):
    """One line, cut on a sentence or a word, no longer than ``limit``."""
    text = " ".join(text.split())
    if len(text) <= limit:
        return text
    window = text[:limit]
    end = window.rfind(". ")
    if end >= 24:
        return window[: end + 1]
    space = window.rfind(" ")
    if space >= 24:
        return window[:space].rstrip(".,;:")
    return window.rstrip()


def _llms():
    global _LLMS
    if _LLMS is None:
        path = ROOT / "scripts" / "gen_llms_txt.py"
        spec = importlib.util.spec_from_file_location("gen_llms_txt", path)
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        _LLMS = module
    return _LLMS


def page_summary(pagename):
    """First sentence of the page source, or the project summary."""
    module = _llms()
    path = module.source_of(pagename)
    if path is not None:
        _title, summary = module.title_and_paragraph(path)
        if summary:
            return clip(summary)
    return clip(module.SUMMARY)


def _at_least(text, title):
    """Keep a real sentence, and make a short one long enough for a description."""
    text = clip(text)
    if len(text) >= 40:
        return text
    extra = "Part of the d-SEAMS manual."
    if text:
        text = text.rstrip(".") + ". " + extra
    else:
        text = "{} is part of the d-SEAMS manual.".format((title or "This page").strip())
    if len(text) < 40:
        text = text.rstrip(".") + ". It documents the engine."
    return clip(text)


def describe_page(pagename, title, doctree):
    """A page-specific description. Shared boilerplate keeps the page title."""
    module = _llms()
    generic = clip(module.SUMMARY)
    specific = ""
    source = page_summary(pagename)
    if source and source != generic and not source.startswith(":") and len(source) >= 12:
        specific = source
    if not specific and doctree is not None:
        for node in _walk_nodes(doctree):
            if _tag(node) != "paragraph":
                continue
            text = clip(node.astext())
            if len(text) >= 12 and not text.startswith(":") and text != generic:
                specific = text
                break
    if specific == generic or not specific:
        specific = "{}. {}".format((title or pagename).strip().rstrip("."), generic)
    return _at_least(specific, title or pagename)


def canonical_url(pagename):
    return "{}/{}.html".format(SITE, pagename)


def markdown_url(pagename):
    return "{}/{}.md".format(SITE, pagename)


def _citation():
    return {
        "@type": "ScholarlyArticle",
        "name": (
            "d-SEAMS: Deferred Structural Elucidation Analysis "
            "for Molecular Simulations"
        ),
        "identifier": DOI_URL,
        "url": DOI_URL,
        "datePublished": "2020-03-20",
        "author": [
            {"@type": "Person", "name": "Rohit Goswami"},
            {"@type": "Person", "name": "Amrita Goswami"},
            {"@type": "Person", "name": "Jayant K. Singh"},
        ],
        "isPartOf": {
            "@type": "PublicationVolume",
            "volumeNumber": "60",
            "isPartOf": {
                "@type": "Periodical",
                "name": "Journal of Chemical Information and Modeling",
            },
        },
    }


def json_ld(pagename, title, description):
    """SoftwareSourceCode, SoftwareApplication, and this WebPage."""
    version = project_version()
    citation = _citation()
    source_languages = ["C++", "C"]
    use_languages = [
        "C++",
        "C",
        "Python",
        "Lua",
        "Fennel",
        "R",
        "Julia",
        "Rust",
        "Fortran",
    ]
    graph = {
        "@context": "https://schema.org",
        "@graph": [
            {
                "@type": "SoftwareSourceCode",
                "name": "d-SEAMS",
                "description": clip(_llms().SUMMARY, 300),
                "version": version,
                "license": LICENSE,
                "codeRepository": REPO,
                "url": SITE + "/",
                "programmingLanguage": source_languages,
                "datePublished": "2026-11-15",
                "citation": citation,
            },
            {
                "@type": "SoftwareApplication",
                "name": "d-SEAMS",
                "applicationCategory": "Science",
                "operatingSystem": "Linux, macOS, Windows",
                "softwareVersion": version,
                "license": LICENSE,
                "codeRepository": REPO,
                "url": SITE + "/",
                "programmingLanguage": use_languages,
                "citation": citation,
                "offers": {"@type": "Offer", "price": "0", "priceCurrency": "USD"},
            },
            {
                "@type": "WebPage",
                "name": title,
                "description": description,
                "url": canonical_url(pagename),
                "isPartOf": SITE + "/",
                "inLanguage": "en",
            },
        ],
    }
    raw = json.dumps(graph, ensure_ascii=True, separators=(",", ":"))
    return raw.replace("<", "\\u003c")


def social_png():
    """1200x630 RGB card. No third-party image library."""
    width, height = 1200, 630
    bg = bytes((12, 92, 86))
    ink = bytes((255, 255, 255))
    band = bytes((183, 228, 216))
    row = bg * width
    canvas = bytearray(row * height)
    for y in range(0, 18):
        start = y * width * 3
        canvas[start : start + width * 3] = band * width
    for y in range(height - 18, height):
        start = y * width * 3
        canvas[start : start + width * 3] = band * width

    def blit(text, scale, y):
        gap = max(scale, 2)
        advance = 5 * scale + gap
        known = [ch for ch in text if ch in _FONT]
        total = len(known) * advance - gap
        x = max(24, (width - total) // 2)
        for ch in known:
            glyph = _FONT[ch]
            for row_i, bits in enumerate(glyph):
                for col in range(5):
                    if not bits & (1 << (4 - col)):
                        continue
                    for dy in range(scale):
                        py = y + row_i * scale + dy
                        if not 0 <= py < height:
                            continue
                        for dx in range(scale):
                            px = x + col * scale + dx
                            if 0 <= px < width:
                                i = (py * width + px) * 3
                                canvas[i : i + 3] = ink
            x += advance

    blit("D-SEAMS", 11, 150)
    blit("NAMES THE ICE IN A TRAJECTORY", 4, 300)
    blit("CHILL+   CAGES   RINGS   QL", 4, 390)
    blit("DOCS.DSEAMS.INFO", 4, 490)
    return _png(width, height, canvas)


def _png(width, height, canvas):
    def chunk(tag, data):
        crc = zlib.crc32(tag + data) & 0xFFFFFFFF
        return struct.pack(">I", len(data)) + tag + data + struct.pack(">I", crc)

    raw = bytearray()
    stride = width * 3
    for y in range(height):
        raw.append(0)
        raw.extend(canvas[y * stride : (y + 1) * stride])
    ihdr = struct.pack(">IIBBBBB", width, height, 8, 2, 0, 0, 0)
    return (
        b"\x89PNG\r\n\x1a\n"
        + chunk(b"IHDR", ihdr)
        + chunk(b"IDAT", zlib.compress(bytes(raw), 9))
        + chunk(b"IEND", b"")
    )


def png_size(data):
    """Return (width, height) from a PNG IHDR, or None."""
    sig = b"\x89PNG\r\n\x1a\n"
    if not data.startswith(sig) or data[12:16] != b"IHDR":
        return None
    width, height = struct.unpack(">II", data[16:24])
    return width, height


def _tag(node):
    return getattr(node, "tagname", node.__class__.__name__)


def _slug(text):
    slug = re.sub(r"[^a-z0-9]+", "-", text.lower()).strip("-")
    return slug or "section"


def _walk_nodes(node):
    yield node
    for child in getattr(node, "children", ()):
        yield from _walk_nodes(child)


def stable_heading_ids(doctree):
    """Give every section an anchor id of letters, digits, ``_`` and ``-``.

    An id docutils already assigned is kept. A missing or odd id is
    replaced, and ``refid`` targets in the same tree follow it.
    """
    used = set()
    renamed = {}
    for node in _walk_nodes(doctree):
        if _tag(node) != "section":
            continue
        title = ""
        for child in node.children:
            if _tag(child) == "title":
                title = child.astext()
                break
        current = list(node.get("ids") or [])
        good = [item for item in current if ID_OK.match(item)]
        chosen = good[0] if good else ""
        if not chosen:
            base = _slug(title)
            chosen = base
            number = 2
            while chosen in used:
                chosen = "{}-{}".format(base, number)
                number += 1
            for old in current:
                if old != chosen:
                    renamed[old] = chosen
            node["ids"] = [chosen]
        used.add(chosen)
    if not renamed:
        return
    for node in _walk_nodes(doctree):
        refid = node.get("refid") if hasattr(node, "get") else None
        if refid in renamed:
            node["refid"] = renamed[refid]


def _escape(text):
    return text.replace("\\", "\\\\").replace("[", "\\[").replace("]", "\\]")


def _md_uri(uri):
    if not uri or uri.startswith(("mailto:", "javascript:")):
        return uri
    path, sep, frag = uri.partition("#")
    if path.endswith(".html"):
        path = path[:-5] + ".md"
    return path + (sep + frag if sep else "")


def _inline(node):
    parts = []
    children = list(getattr(node, "children", ()) or ())
    if not children:
        if _tag(node) in ("#text", "Text") or node.__class__.__name__ == "Text":
            return _escape(node.astext())
        return ""
    for child in children:
        name = child.__class__.__name__
        tag = _tag(child)
        if tag in ("#text", "Text") or name == "Text":
            parts.append(_escape(child.astext()))
        elif tag == "emphasis":
            parts.append("*" + _inline(child) + "*")
        elif tag == "strong":
            parts.append("**" + _inline(child) + "**")
        elif tag == "literal":
            text = child.astext()
            fence = "``" if "`" in text else "`"
            parts.append(fence + text + fence)
        elif tag == "reference":
            label = _inline(child) or _escape(child.astext())
            uri = child.get("refuri") or ""
            if not uri and child.get("refid"):
                uri = "#" + child["refid"]
            uri = _md_uri(uri)
            parts.append("[{}]({})".format(label, uri) if uri else label)
        elif tag == "math":
            parts.append("$" + child.astext() + "$")
        elif tag in ("footnote_reference", "citation_reference"):
            parts.append("[" + child.astext() + "]")
        elif tag == "subscript":
            parts.append("~" + _inline(child) + "~")
        elif tag == "superscript":
            parts.append("^" + _inline(child) + "^")
        elif tag in ("target", "comment", "pending", "problematic", "system_message"):
            continue
        else:
            parts.append(_inline(child))
    return "".join(parts)


def _lang(node):
    lang = node.get("language") or ""
    if lang:
        return str(lang)
    skip = {"code", "code-block", "literal-block", "highlight", "doctest"}
    for name in node.get("classes") or []:
        if name not in skip:
            return name
    return ""


def _fence(node):
    lang = _lang(node)
    body = node.astext().replace("\r\n", "\n")
    return "\n```{}\n{}\n```\n\n".format(lang, body.rstrip("\n"))


def _list(node, level, ordered):
    blocks = []
    index = 1
    for child in node.children:
        if _tag(child) != "list_item":
            continue
        body = _blocks(child, level + 1).strip()
        if not body:
            continue
        prefix = "{}.".format(index) if ordered else "-"
        index += 1
        pad = " " * (len(prefix) + 1)
        lines = body.split("\n")
        joined = prefix + " " + lines[0]
        for line in lines[1:]:
            joined += "\n" + (pad + line if line else "")
        blocks.append(joined)
    if not blocks:
        return ""
    return "\n".join(blocks) + "\n\n"


def _table(node):
    rows = []

    def walk(current):
        if _tag(current) == "row":
            cells = []
            for child in current.children:
                if _tag(child) == "entry":
                    text = " ".join(child.astext().split()).replace("|", "\\|")
                    cells.append(text)
            if cells:
                rows.append(cells)
            return
        for child in getattr(current, "children", ()):
            walk(child)

    walk(node)
    if not rows:
        return ""
    width = max(len(row) for row in rows)
    padded = [row + [""] * (width - len(row)) for row in rows]

    def line(row):
        return "| " + " | ".join(row) + " |"

    out = [line(padded[0]), "| " + " | ".join("---" for _ in range(width)) + " |"]
    out.extend(line(row) for row in padded[1:])
    return "\n".join(out) + "\n\n"


def _admonition(node, level):
    label = _tag(node).replace("_", " ").capitalize()
    body = "".join(_blocks(child, level + 1) for child in node.children).strip()
    lines = ["**{}.** {}".format(label, body.split("\n", 1)[0] if body else "")]
    rest = body.split("\n")[1:] if body else []
    lines.extend(rest)
    return "\n".join("> " + line if line else ">" for line in lines) + "\n\n"


def _blocks(node, level=0):
    tag = _tag(node)
    if tag == "document":
        return "".join(_blocks(child, level) for child in node.children)
    if tag == "section":
        title = None
        rest = []
        for child in node.children:
            if title is None and _tag(child) == "title":
                title = child
            else:
                rest.append(child)
        depth = min(max(level + 1, 1), 6)
        text = _inline(title).strip() if title is not None else "Section"
        ident = ""
        ids = node.get("ids") or []
        if ids and ID_OK.match(ids[0]):
            ident = " {#" + ids[0] + "}"
        body = "".join(_blocks(child, depth) for child in rest)
        return "\n" + ("#" * depth) + " " + text + ident + "\n\n" + body
    if tag in ("paragraph", "compact_paragraph"):
        text = _inline(node).strip()
        return (text + "\n\n") if text else ""
    if tag in ("literal_block", "doctest_block"):
        return _fence(node)
    if tag == "bullet_list":
        return _list(node, level, False)
    if tag == "enumerated_list":
        return _list(node, level, True)
    if tag == "list_item":
        return "".join(_blocks(child, level) for child in node.children)
    if tag == "block_quote":
        body = "".join(_blocks(child, level) for child in node.children).strip()
        if not body:
            return ""
        return "\n".join("> " + line if line else ">" for line in body.split("\n")) + "\n\n"
    if tag in (
        "note",
        "warning",
        "tip",
        "hint",
        "important",
        "caution",
        "attention",
        "danger",
        "error",
        "seealso",
        "admonition",
    ):
        return _admonition(node, level)
    if tag in ("table", "tabular"):
        return _table(node)
    if tag == "image":
        uri = _md_uri(node.get("uri") or "")
        alt = _escape(node.get("alt") or "")
        return "\n![{}]({})\n\n".format(alt, uri)
    if tag == "figure":
        return "".join(_blocks(child, level) for child in node.children)
    if tag == "caption":
        text = _inline(node).strip()
        return ("*" + text + "*\n\n") if text else ""
    if tag == "math_block":
        return "\n$$\n{}\n$$\n\n".format(node.astext().strip())
    if tag == "rubric":
        text = _inline(node).strip()
        return ("**{}**\n\n".format(text)) if text else ""
    if tag == "transition":
        return "\n---\n\n"
    if tag == "desc_signature":
        text = " ".join(node.astext().split())
        return "\n`{}`\n\n".format(text) if text else ""
    if tag in ("definition_list", "definition_list_item", "field_list", "field"):
        return "".join(_blocks(child, level) for child in node.children)
    if tag in ("term", "field_name"):
        text = _inline(node).strip()
        return ("**{}**\n\n".format(text)) if text else ""
    if tag in ("definition", "field_body", "desc", "desc_content", "container", "compound"):
        return "".join(_blocks(child, level) for child in node.children)
    if tag == "line_block":
        lines = [_inline(child).strip() for child in node.children if _tag(child) == "line"]
        return "\n".join(line for line in lines if line) + "\n\n"
    if tag == "raw":
        text = node.astext()
        found = []
        for match in re.finditer(
            r'<div\s+id="(video-[A-Za-z0-9_-]+)"\s*>\s*<h2(?:\s+id="([^"]+)")?>([^<]+)</h2>\s*<p>([^<]+)</p>',
            text,
        ):
            ident = match.group(2) or match.group(1)
            found.append(
                "\n## {} {{#{}}}\n\n{}\n\n".format(
                    match.group(3).strip(),
                    ident,
                    match.group(4).strip(),
                )
            )
        return "".join(found)
    if tag in (
        "comment",
        "target",
        "substitution_definition",
        "pending",
        "system_message",
        "meta",
        "title",
    ):
        return ""
    children = list(getattr(node, "children", ()) or ())
    if children:
        return "".join(_blocks(child, level) for child in children)
    return ""


def page_markdown(doctree, pagename, title):
    url = markdown_url(pagename)
    body = _blocks(doctree).strip() if doctree is not None else ""
    if not body:
        body = "# {}\n".format(title.strip() or pagename)
    text = "<!-- {} -->\n\n{}\n".format(url, body)
    return re.sub(r"\n{3,}", "\n\n", text)


def _git_time(path):
    key = str(path)
    if key in _GIT_TIME:
        return _GIT_TIME[key]
    if path is None or not Path(path).is_file():
        _GIT_TIME[key] = None
        return None
    proc = subprocess.run(
        ["git", "log", "-1", "--format=%cI", "--", key],
        cwd=ROOT,
        check=False,
        capture_output=True,
        text=True,
    )
    found = proc.stdout.strip() or None
    _GIT_TIME[key] = found
    return found


def head_time():
    global _HEAD_TIME
    if _HEAD_TIME is None:
        proc = subprocess.run(
            ["git", "log", "-1", "--format=%cI"],
            cwd=ROOT,
            check=False,
            capture_output=True,
            text=True,
        )
        _HEAD_TIME = proc.stdout.strip() or "2026-11-15T00:00:00+00:00"
    return _HEAD_TIME


def docname_from_loc(loc):
    path = urlsplit(loc).path
    if path.startswith("/"):
        path = path[1:]
    if path.endswith("/"):
        path += "index.html"
    if path.endswith(".html"):
        path = path[:-5]
    return path


def lastmod_for(loc):
    docname = docname_from_loc(loc)
    source = _llms().source_of(docname) if docname else None
    found = _git_time(source) if source is not None else None
    if found:
        return found
    return head_time()


def stamp_sitemap(path, lookup):
    """Add a valid ``lastmod`` to every URL that lacks one."""
    ET.register_namespace("", SITEMAP_NS)
    tree = ET.parse(path)
    root = tree.getroot()
    tag = "{{{}}}".format(SITEMAP_NS)
    for url in list(root):
        if url.tag != tag + "url":
            continue
        loc = None
        last = None
        for child in list(url):
            if child.tag == tag + "loc":
                loc = child
            elif child.tag == tag + "lastmod":
                last = child
        if loc is None or not (loc.text or "").strip():
            continue
        if last is not None and (last.text or "").strip():
            continue
        if last is None:
            last = ET.SubElement(url, tag + "lastmod")
        last.text = lookup(loc.text.strip())
    tree.write(path, xml_declaration=True, encoding="utf-8")


def robots_txt():
    return (
        "User-agent: *\n"
        "Allow: /\n"
        "Content-Signal: search=yes, ai-input=yes, ai-train=yes\n"
        "\n"
        "# Markdown index: {}/llms.txt\n"
        "Sitemap: {}/sitemap.xml\n"
    ).format(SITE, SITE)


def not_found_html():
    version = project_version()
    graph = json_ld("404", "Page not found", "That page is not in the d-SEAMS manual.")
    return """<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <title>Page not found - d-SEAMS</title>
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <meta name="robots" content="noindex">
  <meta name="description" content="That page is not in the d-SEAMS manual.">
  <link rel="canonical" href="{site}/404.html">
  <style>
    body {{ margin: 0; font: 18px/1.5 Georgia, serif; background: #f4fffc; color: #102220; }}
    main {{ max-width: 38rem; margin: 4rem auto; padding: 0 1.5rem; }}
    a {{ color: #0c5c56; }}
    p.mark {{ letter-spacing: 0.08em; text-transform: uppercase; font: 600 13px/1 sans-serif; color: #0c5c56; }}
  </style>
  <script type="application/ld+json">{graph}</script>
</head>
<body>
<main>
  <p class="mark">d-SEAMS {version}</p>
  <h1>Page not found</h1>
  <p>That address is not a page of the manual.</p>
  <ul>
    <li><a href="{site}/index.html">Manual</a></li>
    <li><a href="{site}/quickstart.html">Quickstart</a></li>
    <li><a href="{site}/llms.txt">llms.txt</a></li>
  </ul>
</main>
</body>
</html>
""".format(site=SITE, version=version, graph=graph)


def _captions():
    module = _llms()
    found = {}
    for path in sorted(module.ORG.rglob("*.org")):
        text = path.read_text(encoding="utf-8", errors="replace")
        parent = path.relative_to(module.ORG).with_suffix("").as_posix()
        parent_dir = parent.rsplit("/", 1)[0] if "/" in parent else ""
        title = ""
        for line in text.splitlines():
            if line.upper().startswith("#+TITLE:"):
                title = line.split(":", 1)[1].strip()
                break
        for caption, names in module.toctrees(text):
            label = caption or title
            if not label:
                continue
            for name in names:
                resolved = name.lstrip("/") if name.startswith("/") else "/".join(filter(None, (parent_dir, name)))
                # normalize .. segments
                parts = []
                for part in resolved.split("/"):
                    if part == "..":
                        if parts:
                            parts.pop()
                    elif part and part != ".":
                        parts.append(part)
                found["/".join(parts)] = label
    return found


def _title_map(app):
    titles = {}
    env_titles = getattr(app.env, "titles", {}) or {}
    for name, node in env_titles.items():
        try:
            titles[name] = " ".join(node.astext().split())
        except Exception:
            titles[name] = name
    return titles


def display_title(pagename, titles, captions):
    """Disambiguate pages that share a heading, including inside one section."""
    plain = titles.get(pagename) or pagename
    same = [name for name, text in titles.items() if text == plain]
    if len(same) < 2:
        return plain
    caption = captions.get(pagename) or pagename.split("/", 1)[0]

    def group_of(name):
        return captions.get(name) or name.split("/", 1)[0]

    group = [name for name in same if group_of(name) == caption]
    if len(group) < 2:
        return "{} ({})".format(plain, caption)
    parent = pagename.rsplit("/", 1)[0] if "/" in pagename else ""
    if parent and len([name for name in same if name.startswith(parent + "/")]) == 1:
        kind = parent.rsplit("/", 1)[-1]
        label = {"howto": "How-to", "tutorials": "Tutorial"}.get(kind, kind)
        return "{} ({})".format(plain, label)
    return "{} ({})".format(plain, pagename)


def _on_doctree(app, doctree, docname):
    if getattr(app.builder, "format", None) != "html":
        return
    stable_heading_ids(doctree)


def _on_page(app, pagename, templatename, context, doctree):
    if getattr(app.builder, "format", None) != "html":
        return
    if not hasattr(app, "_dseams_titles"):
        app._dseams_titles = _title_map(app)
        app._dseams_captions = _captions()
    plain = display_title(pagename, app._dseams_titles, app._dseams_captions)
    description = describe_page(pagename, plain, doctree)
    context["title"] = plain
    meta = context.get("meta")
    if not isinstance(meta, dict):
        meta = {}
        context["meta"] = meta
    meta["description"] = description
    meta["canonical_url"] = canonical_url(pagename)
    meta["cover"] = OG_URL
    context["pageurl"] = canonical_url(pagename)
    context["dseams_title"] = "d-SEAMS" if pagename == "index" else plain + " - d-SEAMS"
    context["dseams_description"] = description
    context["dseams_image"] = OG_URL
    context["dseams_markdown"] = markdown_url(pagename)
    context["dseams_jsonld"] = json_ld(pagename, context["dseams_title"], description)
    context["dseams_robots"] = (
        "noindex" if pagename in ("search", "genindex") else "index, follow"
    )
    md_url = markdown_url(pagename)

    def view_markdown(_filename, url=md_url):
        return url

    context["view_source_link"] = view_markdown
    context["page_source_suffix"] = ".md"
    out = Path(app.outdir) / (pagename + ".md")
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(page_markdown(doctree, pagename, plain), encoding="utf-8")


def _llms_full(outdir):
    text = _llms().render()
    parts = ["<!-- Generated from the Markdown twins. Prefer llms.txt, then one page. -->", ""]
    seen = set()
    for match in re.finditer(r"\((https://docs\.dseams\.info/[^)\s]+\.md)\)", text):
        url = match.group(1)
        if url in seen:
            continue
        seen.add(url)
        rel = urlsplit(url).path.lstrip("/")
        path = outdir / rel
        if not path.is_file():
            continue
        parts.append(path.read_text(encoding="utf-8").rstrip())
        parts.append("\n\n---\n")
    body = "\n".join(parts).rstrip() + "\n"
    (outdir / "llms-full.txt").write_text(body, encoding="utf-8")
    (outdir / "llms.txt").write_text(text, encoding="utf-8")


def full_manual_problems(text):
    """Problems with a concatenated manual. Empty when the text is one."""
    problems = []
    if len(text) < 1000:
        problems.append("llms-full.txt is too small to be the manual")
    if "docs.dseams.info/getting-started/quickstart.md" not in text:
        problems.append("llms-full.txt does not contain the quickstart page")
    return problems


def write_manual_for_linkcheck(directory):
    """Write ``llms-full.txt`` for the linkcheck builder.

    The assistants page links ``https://docs.dseams.info/llms-full.txt``.
    The HTML build writes that file from the Markdown twins after it
    renders each page. Linkcheck does not render those twins, and the
    published site does not have the file until this book is deployed,
    so a request to the live URL is a 404. The same writer runs here
    on the source of each page ``llms.txt`` names. Sphinx then checks
    the absolute path as a local file.
    """
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)
    module = _llms()
    seen = set()
    for match in re.finditer(
        r"\((https://docs\.dseams\.info/[^)\s]+\.md)\)", module.render()
    ):
        url = match.group(1)
        if url in seen:
            continue
        seen.add(url)
        rel = urlsplit(url).path.lstrip("/")
        name = rel[:-3] if rel.endswith(".md") else rel
        source = module.source_of(name)
        if source is None:
            continue
        twin = directory / rel
        twin.parent.mkdir(parents=True, exist_ok=True)
        twin.write_text(
            "<!-- {} -->\n\n{}\n".format(url, source.read_text(encoding="utf-8").rstrip()),
            encoding="utf-8",
        )
    _llms_full(directory)
    manual = directory / "llms-full.txt"
    problems = full_manual_problems(manual.read_text(encoding="utf-8"))
    if problems:
        raise RuntimeError("llms-full.txt is not the manual: " + "; ".join(problems))
    return manual


def _on_linkcheck_uri(app, uri):
    """Resolve the one-file manual to the file this build writes."""
    if uri.split("#", 1)[0] != SITE + "/llms-full.txt":
        return None
    builder = getattr(app, "builder", None)
    if builder is not None and getattr(builder, "name", None) != "linkcheck":
        return None
    cached = getattr(_on_linkcheck_uri, "path", None)
    if cached is None:
        manual = write_manual_for_linkcheck(ROOT / "docs" / ".cache" / "llms-linkcheck")
        cached = str(manual.resolve())
        _on_linkcheck_uri.path = cached
    return cached


def _on_finished(app, exception):
    if exception is not None or getattr(app.builder, "format", None) != "html":
        return
    out = Path(app.outdir)
    static = out / "_static"
    static.mkdir(parents=True, exist_ok=True)
    (static / "og-dseams.png").write_bytes(social_png())
    (out / "robots.txt").write_text(robots_txt(), encoding="utf-8")
    (out / "404.html").write_text(not_found_html(), encoding="utf-8")
    _llms_full(out)
    sitemap = out / "sitemap.xml"
    if sitemap.is_file():
        stamp_sitemap(sitemap, lastmod_for)


def connect(app):
    app.connect("doctree-resolved", _on_doctree)
    app.connect("html-page-context", _on_page, priority=800)
    app.connect("build-finished", _on_finished, priority=800)
    app.connect("linkcheck-process-uri", _on_linkcheck_uri)
