#!/usr/bin/env python3
"""Write Sphinx targets for Doxyrest refs that have no anchor.

Doxygen gives some enums and functions a group member id. Doxyrest links
that id from signatures and prose, and the group page does not emit the
anchor. This reads docs/source/xml and writes the missing anchors into
docs/source/api/generated_targets.rst, then lists that file in the API
toctree.
"""

from __future__ import annotations

import re
import sys
import xml.etree.ElementTree as ET
from pathlib import Path


def local(tag: str) -> str:
    return tag.rsplit("}", 1)[-1]


def text_of(el: ET.Element | None) -> str:
    if el is None:
        return ""
    parts = [t.strip() for t in el.itertext() if t and t.strip()]
    return " ".join(parts)


def main() -> int:
    source = Path(sys.argv[1] if len(sys.argv) > 1 else ".").resolve()
    api = source / "api"
    xml_dir = source / "xml"
    if not api.is_dir() or not xml_dir.is_dir():
        print("doxyrest_targets: need api/ and xml/ under", source, file=sys.stderr)
        return 1

    ref_re = re.compile(r":ref:`[^`<]*<doxid-([^>]+)>`")
    explicit_re = re.compile(r"^\.\.\s+_doxid-([^:\s]+):", re.M)
    target_re = re.compile(r":target:`[^`<]*<doxid-([^>]+)>`")

    refs: set[str] = set()
    defined: set[str] = set()
    for path in api.glob("*.rst"):
        if path.name == "generated_targets.rst":
            continue
        text = path.read_text(encoding="utf-8", errors="replace")
        refs.update(ref_re.findall(text))
        defined.update(explicit_re.findall(text))
        defined.update(target_re.findall(text))
    missing = sorted(refs - defined)

    info: dict[str, tuple[str, str]] = {}
    for xml_path in xml_dir.glob("*.xml"):
        try:
            root = ET.parse(xml_path).getroot()
        except ET.ParseError:
            continue
        for el in root.iter():
            if local(el.tag) not in ("memberdef", "compounddef"):
                continue
            ident = el.get("id")
            if not ident or ident not in missing or ident in info:
                continue
            name = ""
            brief = ""
            for child in list(el):
                tag = local(child.tag)
                if tag == "qualifiedname" and child.text:
                    name = child.text.strip()
                elif tag == "name" and child.text and not name:
                    name = child.text.strip()
                elif tag == "title" and child.text and not name:
                    name = child.text.strip()
                elif tag == "briefdescription" and not brief:
                    brief = text_of(child)
            if not name:
                name = ident
            info[ident] = (name, brief)

    lines = [
        "Generated anchors",
        "=================",
        "",
        "Targets for Doxygen ids that the group pages do not anchor.",
        "",
    ]
    for ident in missing:
        name, brief = info.get(ident, (ident, ""))
        lines.append(f".. _doxid-{ident}:")
        lines.append("")
        lines.append(name)
        lines.append("=" * len(name))
        lines.append("")
        if brief:
            lines.append("::")
            lines.append("")
            for row in brief.splitlines() or [brief]:
                lines.append(f"   {row}")
            lines.append("")
    out = api / "generated_targets.rst"
    out.write_text("\n".join(lines).rstrip() + "\n", encoding="utf-8")

    index = api / "index.rst"
    text = index.read_text(encoding="utf-8")
    needle = "generated_targets.rst"
    if needle not in text:
        marker = ".. toctree::\n\t:hidden:\n\n"
        if marker not in text:
            print("doxyrest_targets: API index has no hidden toctree", file=sys.stderr)
            return 1
        text = text.replace(marker, marker + "\tgenerated_targets.rst\n", 1)
        index.write_text(text, encoding="utf-8")
    print(f"doxyrest_targets: {len(missing)} anchors")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
