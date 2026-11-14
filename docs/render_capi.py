#!/usr/bin/env python3
"""Write the C API book from the comments in include/seams.h.

The pages list every SEAMS_API function in this tree. A name that is
only in an older branch is not documented here.

    python3 docs/render_capi.py
"""

import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
HEADER = ROOT / "include" / "seams.h"
OUT = ROOT / "docs" / "source" / "capi"
SYMBOLS = ROOT / "docs" / "capi" / "symbols.txt"

# First matching rule wins. Anything left is "further".
GROUPS = [
    (
        "atoms",
        "Atoms in",
        "Read a frame into an opaque handle.",
        "docs/examples/capi_atoms.c",
        ("seams_system_", "seams_path_"),
    ),
    (
        "neighbours",
        "Neighbours",
        "Build a neighbour list and read it back.",
        "docs/examples/capi_neighbours.c",
        ("seams_neighbours_",),
    ),
    (
        "order",
        "Order parameters",
        "CHILL, CHILL+, and Steinhardt parameters.",
        "docs/examples/capi_order.c",
        ("seams_chill", "seams_ql", "seams_q4q6", "seams_q6bar", "seams_steinhardt"),
    ),
    (
        "rings",
        "Rings",
        "Primitive rings and their histogram.",
        "docs/examples/capi_rings.c",
        ("seams_rings_", "seams_ring_"),
    ),
    (
        "cages",
        "Cages",
        "Hexagonal, cubic, and signature cages.",
        "docs/examples/capi_cages.c",
        ("seams_cages", "seams_signature_cages", "seams_system_symmetry", "seams_prisms", "seams_monolayer"),
    ),
    (
        "errors",
        "Errors",
        "Status codes, the thread-local last error, and version strings.",
        "docs/examples/capi_errors.c",
        (
            "seams_abi",
            "seams_version",
            "seams_status_string",
            "seams_last_error",
            "seams_handle_error",
            "seams_label_name",
            "seams_set_num_threads",
            "seams_get_num_threads",
        ),
    ),
]


def heading(title, char):
    return title + "\n" + (char * len(title)) + "\n"


def rst_text(value):
    return value.replace("\\", "\\\\").replace("*", "\\*").replace("`", "\\`")


def _comment_lines(block):
    lines = []
    for line in block.splitlines():
        line = re.sub(r"^\s*\*\s?", "", line)
        lines.append(line.rstrip())
    return lines


def _parse_comment(block):
    brief = []
    params = []
    returned = []
    prose = []
    current = None
    for line in _comment_lines(block):
        if line.startswith("@brief"):
            current = "brief"
            brief.append(line[len("@brief") :].strip())
            continue
        if line.startswith("@param"):
            current = "param"
            rest = line[len("@param") :].strip()
            name, _, detail = rest.partition(" ")
            params.append([name, detail.strip()])
            continue
        if line.startswith("@return"):
            current = "return"
            returned.append(line[len("@return") :].strip())
            continue
        if line.startswith("@"):
            current = "prose"
            prose.append(line)
            continue
        if line.startswith(" ") and current == "param" and params:
            params[-1][1] = (params[-1][1] + " " + line.strip()).strip()
            continue
        if line.startswith(" ") and current == "brief":
            brief.append(line.strip())
            continue
        if line.startswith(" ") and current == "return":
            returned.append(line.strip())
            continue
        if line.strip() == "":
            current = "prose"
            prose.append("")
            continue
        current = "prose"
        prose.append(line.strip())
    return {
        "brief": " ".join(part for part in brief if part).strip(),
        "params": [(name, text) for name, text in params],
        "returned": " ".join(part for part in returned if part).strip(),
        "prose": "\n".join(prose).strip(),
    }


def parse_header(text):
    found = []
    pattern = re.compile(r"/\*\*(.*?)\*/(.*?)(?=/\*\*|\Z)", re.S)
    # Declarations are not always immediately after the comment: a
    # SEAMS_DEPRECATED line can sit between them. Walk comments and take
    # the next SEAMS_API declaration.
    comments = list(re.finditer(r"/\*\*(.*?)\*/", text, re.S))
    decl = re.compile(
        r"SEAMS_DEPRECATED\(\"([^\"]*)\"\)\s*SEAMS_API\b(.*?);|SEAMS_API\b(.*?);",
        re.S,
    )
    for match in comments:
        tail = text[match.end() : match.end() + 800]
        found_decl = decl.search(tail)
        if found_decl is None or found_decl.start() > 200:
            continue
        if found_decl.group(1) is not None:
            deprecated = found_decl.group(1)
            signature = " ".join(found_decl.group(2).split())
        else:
            deprecated = ""
            signature = " ".join(found_decl.group(3).split())
        name_match = re.search(r"\b(seams_[A-Za-z0-9_]+)\s*\(", signature)
        if name_match is None:
            continue
        record = _parse_comment(match.group(1))
        record["name"] = name_match.group(1)
        record["signature"] = signature
        record["deprecated"] = deprecated
        found.append(record)
    return found


def _group_of(name):
    # Longer, specific prefixes have to win over seams_system_.
    ranked = sorted(GROUPS, key=lambda item: max(len(prefix) for prefix in item[4]), reverse=True)
    for slug, _title, _intro, _example, prefixes in ranked:
        if name.startswith(prefixes):
            return slug
    return "further"


def page(title, intro, example, records, total):
    lines = [heading(title, "="), ""]
    lines.append(intro)
    lines.append("")
    lines.append(
        "Generated from the comments on the ``SEAMS_API`` declarations in "
        "``include/seams.h``. This tree has {} exports. "
        "``meson test docs_capi_symbols`` fails if a declaration is missing "
        "from these pages, and ``meson test docs_capi_examples`` compiles "
        "the example when this page has one.".format(total)
    )
    lines.append("")
    if example:
        lines.append(".. literalinclude:: ../../../" + example)
        lines.append("   :language: c")
        lines.append("")
    for record in records:
        lines.append(".. _capi-{}:".format(record["name"]))
        lines.append("")
        lines.append(heading(record["name"], "-"))
        lines.append("")
        lines.append(".. code-block:: c")
        lines.append("")
        lines.append("   {}".format(record["signature"]))
        lines.append("")
        if record["deprecated"]:
            lines.append("Deprecated. {}.".format(rst_text(record["deprecated"])))
            lines.append("")
        if record["brief"]:
            lines.append(rst_text(record["brief"]))
            lines.append("")
        if record["params"]:
            for name, detail in record["params"]:
                lines.append("``{}``".format(name))
                lines.append("    {}".format(rst_text(detail) if detail else "See the declaration."))
                lines.append("")
        if record["returned"]:
            lines.append("Returns: {}".format(rst_text(record["returned"])))
            lines.append("")
        if record["prose"]:
            for paragraph in record["prose"].split("\n\n"):
                flat = " ".join(paragraph.split())
                if flat:
                    lines.append(rst_text(flat))
                    lines.append("")
    lines.append("")
    return "\n".join(lines)


def main(argv=None):
    del argv
    records = parse_header(HEADER.read_text(encoding="utf-8"))
    names = [record["name"] for record in records]
    if len(names) != len(set(names)):
        raise SystemExit("duplicate SEAMS_API declaration")
    if not names:
        raise SystemExit("no SEAMS_API functions were read from include/seams.h")
    by_group = {slug: [] for slug, *_rest in GROUPS}
    by_group["further"] = []
    for record in records:
        by_group[_group_of(record["name"])].append(record)
    OUT.mkdir(parents=True, exist_ok=True)
    SYMBOLS.parent.mkdir(parents=True, exist_ok=True)
    ordered = []
    index = [heading("C API", "="), ""]
    index.append(
        "Every ``SEAMS_API`` function in ``include/seams.h``. "
        "The list is generated from that header, so a function that is not "
        "declared there is not documented here. "
        "``meson test docs_capi_symbols`` reads the header again and fails "
        "if one of those declarations is absent from these pages."
    )
    index.append("")
    index.append(".. toctree::")
    index.append("   :maxdepth: 1")
    index.append("")
    blocks = list(GROUPS) + [
        (
            "further",
            "Further analyses",
            "Brine fronts, topology, and the remaining exports.",
            None,
            (),
        )
    ]
    for slug, title, intro, example, _prefixes in blocks:
        group = by_group[slug]
        if not group:
            continue
        ordered.extend(record["name"] for record in group)
        (OUT / (slug + ".rst")).write_text(
            page(title, intro, example, group, len(names)), encoding="utf-8"
        )
        index.append("   " + slug)
    index.append("   merge")
    index.append("")
    (OUT / "index.rst").write_text("\n".join(index) + "\n", encoding="utf-8")
    SYMBOLS.write_text("\n".join(ordered) + "\n", encoding="utf-8")
    merge = "\n".join(
        [
            heading("Where the API lives", "="),
            "",
            "``include/seams.h`` in this tree is the C API. "
            "``docs/render_capi.py`` reads that header and writes these pages "
            "and ``docs/capi/symbols.txt``.",
            "",
            "Regenerate them with:",
            "",
            ".. code-block:: bash",
            "",
            "   python3 docs/render_capi.py",
            "",
            "``meson test docs_capi_symbols`` fails if the header and the pages disagree. "
            "``meson test docs_capi_examples`` compiles and runs the six examples. "
            "Those examples call only functions the header declares.",
            "",
        ]
    )
    (OUT / "merge.rst").write_text(merge, encoding="utf-8")
    print("wrote", len(ordered), "functions from", HEADER)
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
