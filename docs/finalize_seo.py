#!/usr/bin/env python3
"""Complete metadata for non-Sphinx HTML and write a whole-site sitemap."""

from __future__ import annotations

import html
import json
import sys
import tempfile
import xml.etree.ElementTree as ET
from pathlib import Path

from check_seo import BASE_URL, MetadataParser

DESCRIPTION = (
    "Documentation for d-SEAMS, a molecular simulation trajectory analysis engine "
    "for ice, clathrates, rings, cages, and structural order."
)


def canonical_url(path: Path, root: Path) -> str:
    return BASE_URL + path.relative_to(root).as_posix()


def _json_ld() -> str:
    data = {
        "@context": "https://schema.org",
        "@graph": [
            {
                "@type": "SoftwareApplication",
                "@id": BASE_URL + "#application",
                "name": "d-SEAMS",
                "description": DESCRIPTION,
                "applicationCategory": "ScientificApplication",
                "operatingSystem": "Linux, macOS",
                "url": BASE_URL,
                "softwareVersion": "2.14.0",
                "license": "https://spdx.org/licenses/GPL-3.0-or-later.html",
                "citation": "https://doi.org/10.1021/acs.jcim.0c00031",
            },
            {
                "@type": "SoftwareSourceCode",
                "@id": "https://github.com/d-SEAMS/seams-core#source",
                "name": "seams-core",
                "codeRepository": "https://github.com/d-SEAMS/seams-core",
                "programmingLanguage": ["C++", "C", "Python"],
                "license": "https://spdx.org/licenses/GPL-3.0-or-later.html",
                "citation": "https://doi.org/10.1021/acs.jcim.0c00031",
            },
        ],
    }
    return json.dumps(data, ensure_ascii=True, separators=(",", ":")).replace("<", "\\u003c")


def complete_page(path: Path, root: Path) -> str:
    source = path.read_text(encoding="utf-8", errors="replace")
    parser = MetadataParser()
    parser.feed(source)
    title = "".join(parser.title).strip()
    if not title:
        raise ValueError(f"{path.relative_to(root)} has no title")
    url = canonical_url(path, root)
    additions = []
    if not any(parser.descriptions):
        additions.append(f'<meta name="description" content="{html.escape(DESCRIPTION, quote=True)}">')
    if not parser.canonicals:
        additions.append(f'<link rel="canonical" href="{html.escape(url, quote=True)}">')
    og_values = {
        "og:title": title,
        "og:description": DESCRIPTION,
        "og:image": "https://dseams.info/images/og-card.png",
        "og:url": url,
        "og:type": "website",
    }
    for prop, value in og_values.items():
        if not any(parser.open_graph.get(prop, [])):
            additions.append(
                f'<meta property="{prop}" content="{html.escape(value, quote=True)}">'
            )
    if 'name="twitter:card"' not in source:
        additions.append('<meta name="twitter:card" content="summary">')
    if 'name="twitter:title"' not in source:
        additions.append(
            f'<meta name="twitter:title" content="{html.escape(title, quote=True)}">'
        )
    if 'name="twitter:description"' not in source:
        additions.append(
            f'<meta name="twitter:description" content="{html.escape(DESCRIPTION, quote=True)}">'
        )
    if 'name="twitter:image"' not in source:
        additions.append(
            '<meta name="twitter:image" '
            'content="https://dseams.info/images/og-card.png">'
        )
    if 'type="application/ld+json"' not in source:
        additions.append(f'<script type="application/ld+json">{_json_ld()}</script>')
    if additions:
        marker = source.lower().find("</head>")
        if marker < 0:
            raise ValueError(f"{path.relative_to(root)} has no closing head element")
        source = source[:marker] + "\n".join(additions) + "\n" + source[marker:]
        path.write_text(source, encoding="utf-8")
    return url


def finalize(root: Path) -> None:
    pages = sorted(root.rglob("*.html")) if root.is_dir() else []
    if not pages:
        raise ValueError(f"{root}: no built HTML pages found")
    urls = [complete_page(page, root) for page in pages]
    ET.register_namespace("", "http://www.sitemaps.org/schemas/sitemap/0.9")
    urlset = ET.Element("{http://www.sitemaps.org/schemas/sitemap/0.9}urlset")
    for url in urls:
        node = ET.SubElement(urlset, "{http://www.sitemaps.org/schemas/sitemap/0.9}url")
        ET.SubElement(node, "{http://www.sitemaps.org/schemas/sitemap/0.9}loc").text = url
    ET.ElementTree(urlset).write(
        root / "sitemap.xml", encoding="utf-8", xml_declaration=True
    )


def self_test() -> None:
    with tempfile.TemporaryDirectory() as tmp:
        root = Path(tmp)
        page = root / "api-cpp" / "index.html"
        page.parent.mkdir()
        page.write_text("<html><head><title>C++ API</title></head></html>", encoding="utf-8")
        finalize(root)
        source = page.read_text(encoding="utf-8")
        required = (
            "canonical",
            "og:title",
            "og:image",
            "twitter:card",
            "twitter:image",
            "application/ld+json",
        )
        if not all(item in source for item in required):
            raise ValueError("self-test did not add all metadata")
        tree = ET.parse(root / "sitemap.xml")
        if not any(node.text == canonical_url(page, root) for node in tree.iter()):
            raise ValueError("self-test sitemap omitted the page")


def main(argv: list[str]) -> int:
    try:
        if argv[1:] == ["--self-test"]:
            self_test()
        elif len(argv) == 2:
            finalize(Path(argv[1]))
        else:
            print("usage: docs/finalize_seo.py BUILD_DIR | --self-test", file=sys.stderr)
            return 2
    except (OSError, ValueError, ET.ParseError) as exc:
        print(exc, file=sys.stderr)
        return 1
    print("SEO metadata finalized")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
