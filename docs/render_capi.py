#!/usr/bin/env python3
"""Write the C API book from Doxygen XML of include/seams.h.

The XML is produced from capi/polish-2.14 at 0127f41e, not from the
docs branch. Pass the xml/ directory:

    python3 docs/render_capi.py /path/to/doxygen/xml
"""

import sys
import xml.etree.ElementTree as ET
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "docs" / "source" / "capi"
SYMBOLS = ROOT / "docs" / "capi" / "symbols.txt"

# Task groups. Every SEAMS_API function is listed once.
GROUPS = [
    (
        "atoms",
        "Atoms in",
        "Read a frame into an opaque handle.",
        "docs/examples/capi_atoms.c",
        [
            "seams_system_create",
            "seams_system_destroy",
            "seams_system_set_atoms",
            "seams_system_set_frame",
            "seams_system_set_ortho",
            "seams_system_read",
            "seams_system_count",
            "seams_system_natoms",
            "seams_path_nframes",
        ],
    ),
    (
        "neighbours",
        "Neighbours",
        "Build a neighbour list and read it back.",
        "docs/examples/capi_neighbours.c",
        [
            "seams_neighbours_build",
            "seams_neighbours_destroy",
            "seams_neighbours_count",
            "seams_neighbours_nrows",
            "seams_neighbours_csr",
            "seams_neighbours_copy",
            "seams_neighbours_knn",
        ],
    ),
    (
        "order",
        "Order parameters",
        "CHILL, CHILL+, and Steinhardt parameters.",
        "docs/examples/capi_order.c",
        [
            "seams_chill",
            "seams_chill_plus",
            "seams_chill_plus_n",
            "seams_ql",
            "seams_q4q6",
            "seams_q4q6_n",
            "seams_q6bar",
            "seams_steinhardt",
        ],
    ),
    (
        "rings",
        "Rings",
        "Primitive rings and their histogram.",
        "docs/examples/capi_rings.c",
        [
            "seams_rings_find",
            "seams_rings_destroy",
            "seams_rings_count",
            "seams_rings_nrings",
            "seams_rings_histogram",
            "seams_rings_histogram_n",
            "seams_rings_labels",
            "seams_rings_labels_n",
            "seams_rings_vertices",
            "seams_rings_fraction",
            "seams_ring_symmetry_name",
        ],
    ),
    (
        "cages",
        "Cages",
        "Hexagonal and cubic cage membership.",
        "docs/examples/capi_cages.c",
        [
            "seams_cages_find",
            "seams_cages",
            "seams_cages_destroy",
            "seams_cages_counts",
            "seams_cages_labels",
            "seams_cages_signature",
            "seams_signature_cages",
            "seams_system_symmetry",
            "seams_prisms",
            "seams_monolayer",
        ],
    ),
    (
        "errors",
        "Errors",
        "Status codes, the thread-local last error, and version strings.",
        "docs/examples/capi_errors.c",
        [
            "seams_abi",
            "seams_version",
            "seams_version_major",
            "seams_version_minor",
            "seams_version_patch",
            "seams_version_string",
            "seams_status_string",
            "seams_last_error",
            "seams_last_error_message",
            "seams_handle_error",
            "seams_label_name",
            "seams_set_num_threads",
            "seams_get_num_threads",
        ],
    ),
    (
        "further",
        "Further analyses",
        "Brine fronts, topology, and the remaining exports.",
        None,
        [
            "seams_brine_analyse",
            "seams_brine_destroy",
            "seams_brine_fronts",
            "seams_brine_pockets",
            "seams_brine_profile",
            "seams_brine_rejection",
            "seams_contact_pairs",
            "seams_coordination",
            "seams_density_axis",
            "seams_domains",
            "seams_fingerprint",
            "seams_fingerprint_atom",
            "seams_fingerprint_destroy",
            "seams_fingerprint_key",
            "seams_fingerprint_rings",
            "seams_hbonds",
            "seams_rdf",
            "seams_running_cn",
            "seams_site_table",
        ],
    ),
]


def plain(node):
    if node is None:
        return ""
    parts = [node.text or ""]
    for child in node:
        if child.tag in ("para", "parameteritem", "simplesect", "listitem", "parameterdescription"):
            parts.append(" ")
        parts.append(plain(child))
        if child.tag == "parametername":
            parts.append(" ")
        parts.append(child.tail or "")
    return "".join(parts)


def load_functions(xml_dir):
    found = {}
    for path in sorted(Path(xml_dir).glob("*.xml")):
        try:
            root = ET.parse(path).getroot()
        except ET.ParseError:
            continue
        for member in root.iter("memberdef"):
            if member.get("kind") != "function":
                continue
            name = (member.findtext("name") or "").strip()
            if not name.startswith("seams_"):
                continue
            definition = (member.findtext("definition") or name).strip()
            args = (member.findtext("argsstring") or "").strip()
            brief = " ".join(plain(member.find("briefdescription")).split())
            detail = " ".join(plain(member.find("detaileddescription")).split())
            params = []
            for param in member.findall("param"):
                decl = (param.findtext("declname") or "").strip()
                if decl:
                    params.append(decl)
            current = found.get(name)
            record = {
                "name": name,
                "signature": (definition + args).replace("  ", " "),
                "brief": brief,
                "detail": detail,
                "params": params,
            }
            if current is None or len(record["detail"]) > len(current["detail"]):
                found[name] = record
    return found


def rst_text(value):
    """Keep Doxygen prose from opening an RST emphasis or interpreted text."""
    return value.replace("\\", "\\\\").replace("*", "\\*").replace("`", "\\`")


def heading(title, char):
    return title + "\n" + (char * len(title)) + "\n"


def page(slug, title, intro, example, names, functions):
    lines = [heading(title, "="), ""]
    lines.append(intro)
    lines.append("")
    lines.append(
        "Generated from the Doxygen XML of ``include/seams.h`` at "
        "``capi/polish-2.14`` commit ``0127f41efc3f44cefdcf532bdea3bba36b05176d``. "
        "The C sources of that branch are not part of this docs tree."
    )
    lines.append("")
    if example:
        lines.append("``meson test docs_capi_examples`` compiles and runs this example.")
        lines.append("")
        lines.append(".. literalinclude:: ../../../" + example)
        lines.append("   :language: c")
        lines.append("")
    for name in names:
        record = functions[name]
        lines.append(".. _capi-{}:".format(name))
        lines.append("")
        lines.append(heading(name, "-"))
        lines.append("")
        lines.append(".. code-block:: c")
        lines.append("")
        lines.append("   {}".format(record["signature"]))
        lines.append("")
        if record["brief"]:
            lines.append(rst_text(record["brief"]))
            lines.append("")
        if record["detail"] and record["detail"] != record["brief"]:
            lines.append(rst_text(record["detail"]))
            lines.append("")
    lines.append("")
    return "\n".join(lines)


def main(argv):
    if len(argv) != 1:
        print("usage: render_capi.py XML_DIR", file=sys.stderr)
        return 1
    functions = load_functions(argv[0])
    wanted = [name for _slug, _title, _intro, _example, names in GROUPS for name in names]
    missing = [name for name in wanted if name not in functions]
    extra = sorted(set(functions) - set(wanted))
    if missing or extra or len(wanted) != len(set(wanted)):
        print("missing", missing, file=sys.stderr)
        print("extra", extra, file=sys.stderr)
        return 1
    OUT.mkdir(parents=True, exist_ok=True)
    SYMBOLS.parent.mkdir(parents=True, exist_ok=True)
    SYMBOLS.write_text("\n".join(wanted) + "\n", encoding="utf-8")
    index = [heading("C API", "="), ""]
    index.append(
        "Every function exported by ``src/libseams.map`` (the ``seams_*`` "
        "version script) on ``capi/polish-2.14``. The concrete list is the "
        "77 ``SEAMS_API`` functions in ``include/seams.h`` at "
        "``0127f41e``. ``meson test docs_capi_symbols`` fails if one of "
        "them is absent from these pages."
    )
    index.append("")
    index.append(".. toctree::")
    index.append("   :maxdepth: 1")
    index.append("")
    for slug, title, intro, example, names in GROUPS:
        (OUT / (slug + ".rst")).write_text(
            page(slug, title, intro, example, names, functions), encoding="utf-8"
        )
        index.append("   " + slug)
    index.append("   merge")
    index.append("")
    (OUT / "index.rst").write_text("\n".join(index) + "\n", encoding="utf-8")
    print("wrote", len(wanted), "functions")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
