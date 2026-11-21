#!/usr/bin/env python3
# Copyright (c) d-SEAMS developers
"""Fill missing contract lines on the C and C++ API reference pages.

Existing purpose, parameter, unit, and example lines stay. On the C page,
sections rendered from include/seams.h by docs/capi_docs.py keep their
header fields and gain any missing Ownership, Errors, or Example lines,
while non-function C symbols follow include/seams.h. On the C++ language
ABI page, existing curated section order and Return, Ownership, and Error
lines stay, and any missing or placeholder fields and new symbols are
filled from bindings/cxx/seams_lang.h. Nothing outside the symbol block
is removed.

Only the C and C++ pages are written here. The binding repositories keep
their own language references.
"""

import argparse
import difflib
import re
import sys
from pathlib import Path

import check_c_symbols

ROOT = Path(__file__).resolve().parents[1]
API = ROOT / "docs" / "api"
BEGIN = "<!-- API_SYMBOLS_BEGIN -->"
END = "<!-- API_SYMBOLS_END -->"
PAGES = {
    "c": "c.md",
    "cxx": "cpp.md",
}
EXAMPLES = {
    "c": "docs/examples/c_api_tour.c",
    "cxx": "docs/examples/minimal_cpp.cpp",
}
STATUS = (
    "SEAMS_ERR_NULL, SEAMS_ERR_RANGE, SEAMS_ERR_BOX, SEAMS_ERR_ARG, "
    "SEAMS_ERR_STATE, SEAMS_ERR_NOMEM, or SEAMS_ERR_INTERNAL. "
    "seams_last_error() holds the detail, and seams_handle_error(handle) "
    "when the call took a handle."
)


def _clean(raw):
    lines = []
    for line in raw.splitlines():
        line = re.sub(r"^\s*\*\s?", "", line).strip()
        if line:
            lines.append(line)
    return " ".join(lines)


def _first_sentence(text):
    text = " ".join(text.split())
    match = re.match(r"(.+?\.)(?:\s|$)", text)
    return (match.group(1) if match else text).strip()


def _slice(text, label):
    match = re.search(
        label + r"\s+(.*?)(?=\s+@|\s+Thread safety:|\s+Ownership:|$)",
        text,
    )
    if match is None:
        return ""
    return _first_sentence(match.group(1)).rstrip(".")


def _args(blob):
    if blob.strip() in ("", "void"):
        return []
    parts = []
    depth = 0
    buf = []
    for char in blob:
        if char == "(":
            depth += 1
        elif char == ")":
            depth -= 1
        if char == "," and depth == 0:
            parts.append("".join(buf))
            buf = []
        else:
            buf.append(char)
    parts.append("".join(buf))
    names = []
    skip = {"const", "int", "double", "void", "char", "size_t", "unsigned", "struct"}
    for part in parts:
        part = part.replace("[9]", "").replace("[]", "").strip()
        match = re.search(r"([A-Za-z_][A-Za-z0-9_]*)\s*$", part)
        if match and match.group(1) not in skip:
            names.append(match.group(1))
    return names


def _sentence(text):
    text = " ".join(text.split())
    if not text:
        return ""
    if text[-1] not in ".!?":
        text += "."
    return text


def _c_contract(name, ret, args, comment):
    ret = " ".join(ret.split())
    brief = _slice(comment, "@brief") or _first_sentence(comment).rstrip(".")
    params = _args(args)
    param_text = "; ".join(params) if params else "none"
    returned = _slice(comment, "@return")
    owned = ""
    own = re.search(r"Ownership:\s+(.*?)(?=\s+Thread safety:|$)", comment)
    if own:
        owned = _first_sentence(own.group(1))
    elif "owned by" in comment:
        for part in re.split(r"(?<=\.)\s+", comment):
            if "owned by" in part:
                owned = _first_sentence(part)
                break
    out = re.findall(r"(seams_[A-Za-z0-9_]+)\s*\*\*", args)
    if ret == "void" and (name.endswith("_destroy") or name.endswith("_free")):
        returns = "Nothing."
        ownership = "Frees the caller-owned handle. delete of a null pointer is a no-op."
        errors = "None."
        units = "handle"
    elif ret.startswith("const char"):
        returns = _sentence(returned or "Pointer to library text")
        ownership = "The library owns the string. The caller does not free it."
        errors = "None."
        units = "text"
    elif ret == "seams_status":
        returns = _sentence(returned or "seams_status. SEAMS_OK (0) is success")
        if not owned:
            if out:
                ownership = (
                    "On success the caller owns the new "
                    + out[0]
                    + " and frees it with the matching destroy function. "
                    "Buffers the caller passed stay with the caller."
                )
            elif "NULL" in comment and "query" in comment:
                ownership = (
                    "The caller owns every buffer. A NULL buffer with a length "
                    "pointer queries the required count and returns SEAMS_OK. "
                    "A short buffer returns SEAMS_ERR_RANGE and writes nothing."
                )
            else:
                ownership = "The caller owns every buffer it passes. This call does not allocate a handle the caller must free."
        else:
            ownership = _sentence(owned)
        errors = STATUS
        units = "cutoff and lengths in angstrom" if "cutoff" in args or "rmax" in args else "as named by the parameters"
    elif "*" in ret and ret.startswith("seams_"):
        returns = "A new handle, or NULL when allocation fails."
        ownership = "The caller owns the handle and frees it with the matching destroy function."
        errors = "NULL. seams_last_error() then reads 'out of memory'."
        units = "handle"
    elif ret in {"int", "size_t"}:
        returns = _sentence(returned or "The integer result")
        ownership = "No allocation."
        errors = "None."
        units = "dimensionless"
    else:
        returns = _sentence(returned or ret)
        ownership = _sentence(owned or "No allocation.")
        errors = "None."
        units = "as named by the parameters"
    return {
        "purpose": _sentence(brief) + " Declared in include/seams.h.",
        "parameters": param_text,
        "returns": returns,
        "ownership": ownership,
        "errors": errors,
        "units": units,
    }


def c_records():
    text = (ROOT / "include" / "seams.h").read_text()
    found = {}
    for match in re.finditer(r"/\*\*(.*?)\*/", text, re.S):
        comment = _clean(match.group(1))
        after = text[match.end() :]
        nxt = re.search(r"/\*\*(?![*<])", after)
        if nxt:
            after = after[: nxt.start()]
        after = re.sub(r"SEAMS_DEPRECATED\([^)]*\)\s*", "", after)
        func = re.match(
            r"\s*SEAMS_API\s+(.*?)(seams_[A-Za-z0-9_]+)\s*\(([^;]*)\)\s*;",
            after,
            re.S,
        )
        if func:
            found[func.group(2)] = _c_contract(func.group(2), func.group(1), func.group(3), comment)
            continue
        typed = re.match(r"\s*typedef struct (\w+) (\w+)\s*;", after)
        if typed:
            name = typed.group(2)
            found[name] = {
                "purpose": _sentence(comment) + " Declared in include/seams.h.",
                "parameters": "none",
                "returns": "The opaque type.",
                "ownership": "The caller creates the object and frees it with the matching destroy function.",
                "errors": "None. This is a type.",
                "units": "handle",
            }
            continue
        enum = re.match(r"\s*typedef enum (\w+) \{(.*?)\}\s*(\w+)\s*;", after, re.S)
        if enum:
            kind = enum.group(3)
            found[kind] = {
                "purpose": _sentence(comment or kind) + " Declared in include/seams.h.",
                "parameters": "none",
                "returns": "The enumeration type.",
                "ownership": "No allocation.",
                "errors": "None. This is a type.",
                "units": "dimensionless",
            }
            for item, raw, note in re.findall(
                r"(SEAMS_[A-Z0-9_]+)\s*=\s*([^,/\n]+),?\s*(?:/\*\*<\s*(.*?)\*/)?",
                enum.group(2),
            ):
                purpose = _sentence((note or comment or item).strip())
                found[item] = {
                    "purpose": purpose + " Declared in include/seams.h.",
                    "parameters": "none",
                    "returns": "The enumerator value " + " ".join(raw.split()) + ".",
                    "ownership": "No allocation.",
                    "errors": "None.",
                    "units": "dimensionless code",
                }
    for name, purpose, value in (
        ("SEAMS_ABI", "ABI number compiled into this header.", "1"),
        ("SEAMS_DEPRECATED", "Mark a shipped declaration and name its replacement.", "attribute"),
        ("SEAMS_VERSION_MAJOR", "Major component of this release.", "integer"),
        ("SEAMS_VERSION_MINOR", "Minor component of this release.", "integer"),
        ("SEAMS_VERSION_PATCH", "Patch component of this release.", "integer"),
        ("SEAMS_VERSION_STRING", "Release text written by the build.", "text"),
        (
            "SEAMS_VERSION_NUMBER",
            "Packed version, (major << 16) | (minor << 8) | patch.",
            "integer",
        ),
    ):
        found.setdefault(
            name,
            {
                "purpose": purpose + " Declared in include/seams.h.",
                "parameters": "message" if name == "SEAMS_DEPRECATED" else "none",
                "returns": "The macro value (" + value + ").",
                "ownership": "No allocation.",
                "errors": "None.",
                "units": "text" if value == "text" else "dimensionless",
            },
        )
    found["seams_lammps_bounds"] = {
        "purpose": "Fill ox, oy, oz, lx, ly, lz, xy, xz, yz from LAMMPS box bounds. Declared in include/seams_lammps_box.h.",
        "parameters": "lo; hi; tilt; box9",
        "returns": "Nothing. The nine cell numbers are written to box9.",
        "ownership": "The caller owns lo, hi, tilt, and box9. The function does not allocate.",
        "errors": "None. The caller supplies three-vectors and a nine-number buffer.",
        "units": "lengths in angstrom",
    }
    return found


LUA_PURPOSE = {
    "system": "Build a cloud from coordinates and a cell, or an empty cloud when coordinates are omitted.",
    "set_num_threads": "Set the OpenMP width for the next call.",
    "num_threads": "Report the OpenMP width the next call will use.",
    "version": "Return the library version text.",
    "label_name": "Map a CHILL code to its name.",
    "read": "Read one frame from a trajectory path.",
    "nframes": "Count frames in a trajectory path.",
    "neighbors": "Build a cutoff, mutual, or seeded neighbour list.",
    "knn": "Build a k-nearest neighbour list.",
    "chill_plus": "Write one CHILL+ label per atom.",
    "chill": "Write one original CHILL label per atom.",
    "q4q6": "Write local q4 and q6 on the cutoff shell.",
    "steinhardtQl": "Write local Steinhardt q_l, and the neighbour average when asked.",
    "rings": "Return primitive-ring sizes, symmetry, a histogram, and the six-ring fraction.",
    "ringNetwork": "Return the same ring census with an explicit maximum size.",
    "getPrimitiveRings": "Return the atom indices of one primitive ring.",
    "cages": "Score hexagonal and double-diamond cage membership.",
    "cages_by_signature": "Count closed cages whose ring census matches a spec.",
    "rdf": "Compute a partial radial distribution.",
    "running_cn": "Integrate the partial radial distribution into a running coordination number.",
    "cn": "Return the coordination number at the requested cutoff.",
    "hbonds": "Count hydrogen bonds on one frame.",
    "density": "Return the number density along one cell axis.",
    "fingerprint": "Build the topological fingerprint of the bonded graph.",
    "site_table": "Assign one site kind per atom from a site-table spec.",
    "pairs": "Count mutual nearest unlike contacts of a site table.",
    "domain": "Measure the largest polar or apolar domain.",
    "prisms": "Count perfect and imperfect prisms by ring size.",
    "monolayer": "Histogram primitive rings in a monolayer.",
    "core": "The C module table the Lua functions call.",
    "cloud.close": "Free the system stored on a cloud.",
    "cloud.box": "Return the nine LAMMPS cell numbers stored on a cloud.",
    "cloud.nop": "Atom count stored on a cloud when it was built.",
}


def _cxx_comment(raw_blocks, name):
    parts = []
    for raw in raw_blocks:
        if "@file" in raw:
            continue
        cleaned = _clean(raw)
        cleaned = re.sub(r"@ingroup\s+\S+", "", cleaned)
        cleaned = re.sub(r"@brief\s+", "", cleaned)
        cleaned = " ".join(cleaned.split())
        if cleaned and cleaned.rstrip(".") != name:
            parts.append(cleaned)
    text = " ".join(parts)
    return _first_sentence(text) if text else ""


def _cxx_contract(name, ret, args, comment):
    ret = " ".join(ret.split())
    params = _args(args)
    param_text = "; ".join(params) if params else "none"
    if ret == "void" and (name.endswith("_free") or name.endswith("_destroy")):
        returns = "Nothing."
        ownership = "Frees the caller-owned frame. A null pointer is a no-op."
        errors = "None."
        units = "handle"
    elif ret == "void":
        returns = "Nothing."
        ownership = "The caller owns every buffer it passed. The frame stays with its owner."
        errors = "None."
        units = "as named by the parameters"
    elif ret in ("seams_frame *", "seams_frame*"):
        returns = "A new seams_frame, or NULL when the call cannot build one."
        ownership = "The caller owns the frame and frees it with seams_frame_free."
        errors = "NULL. seams_last_error() names the reason."
        units = "coordinates and box lengths in angstrom"
    elif ret.startswith("const char"):
        returns = "A pointer to library text."
        ownership = "The library owns the string. The caller does not free it."
        errors = "None."
        units = "text"
    elif ret.startswith("const ") and "*" in ret:
        returns = "A pointer into the frame."
        ownership = (
            "The frame owns the memory. It stays valid until seams_frame_free "
            "or a call that replaces the coordinates."
        )
        errors = "NULL when the frame is null."
        units = "as named by the parameters"
    elif ret == "int" and "seams_frame" in args:
        returns = "Zero on success. A non-zero value is the seams_status code."
        ownership = "The caller owns every buffer it passed. The frame stays with its owner."
        errors = "A non-zero code. seams_last_error() holds the detail."
        units = "cutoff in angstrom" if "cutoff" in args else "as named by the parameters"
    elif ret == "int" and "*" in args:
        returns = "Zero on success. A non-zero value reports failure."
        ownership = "The caller owns every input and output buffer."
        errors = "A non-zero code. seams_last_error() holds the detail."
        units = "dimensionless indices and counts"
    else:
        returns = "The integer result." if ret == "int" else _sentence(ret)
        ownership = "No allocation."
        errors = "None."
        units = "dimensionless"
    purpose = (
        _sentence(comment) + " Declared in bindings/cxx/seams_lang.h."
        if comment
        else "Declared in bindings/cxx/seams_lang.h."
    )
    return {
        "purpose": purpose,
        "parameters": param_text,
        "returns": returns,
        "ownership": ownership,
        "errors": errors,
        "units": units,
    }


def cxx_records():
    text = (ROOT / "bindings" / "cxx" / "seams_lang.h").read_text()
    text = re.sub(
        r"#ifndef\s+DOXYGEN_SHOULD_SKIP_THIS.*?#endif",
        "",
        text,
        flags=re.S,
    )
    found = {}
    token = re.compile(
        r"(?P<doc>/\*\*(?:(?!\*/).)*?\*/)"
        r"|(?P<opaque>typedef\s+struct\s+(?P<oname>\w+)\s+(?P<oalias>\w+)\s*;)"
        r"|(?P<struct>typedef\s+struct\s+(?P<sname>\w+)\s*\{(?P<sbody>[^}]*)\}\s*(?P<salias>\w+)\s*;)"
        r"|(?P<func>(?P<fret>const\s+\w+\s*\*|seams_frame\s*\*|void|int)\s*(?P<fname>seams_[A-Za-z0-9_]+)\s*\((?P<fargs>[^;]*?)\)\s*;)",
        re.S,
    )
    pending = []
    last_end = 0
    for match in token.finditer(text):
        gap = text[last_end : match.start()]
        if gap.strip():
            pending = []
        last_end = match.end()
        kind = match.lastgroup
        if kind == "doc":
            pending.append(match.group("doc")[3:-2])
            continue
        blocks = pending
        pending = []
        if kind == "opaque":
            name = match.group("oalias")
            note = _cxx_comment(blocks, name)
            found[name] = {
                "purpose": (_sentence(note) + " " if note else "") + "Declared in bindings/cxx/seams_lang.h.",
                "parameters": "none",
                "returns": "The type itself.",
                "ownership": "No allocation. This is a type.",
                "errors": "None.",
                "units": "handle",
            }
        elif kind == "struct":
            name = match.group("salias")
            note = _cxx_comment(blocks, name)
            fields = [
                fmatch.group(1)
                for line in match.group("sbody").splitlines()
                for fmatch in [re.search(r"([A-Za-z_][A-Za-z0-9_]*)\s*;\s*$", line)]
                if fmatch
            ]
            found[name] = {
                "purpose": (_sentence(note) + " " if note else "") + "Declared in bindings/cxx/seams_lang.h.",
                "parameters": "; ".join(fields) if fields else "none",
                "returns": "The type itself.",
                "ownership": "No allocation. This is a type.",
                "errors": "None.",
                "units": "as named by the parameters",
            }
        elif kind == "func":
            name = match.group("fname")
            if name in found:
                continue
            found[name] = _cxx_contract(
                name,
                match.group("fret"),
                match.group("fargs"),
                _cxx_comment(blocks, name),
            )
    return found


def semantics(lang, name, purpose, parameters):
    """Return, ownership, and error lines for a symbol that already has a purpose."""
    if lang == "c":
        record = c_records().get(name)
        if record:
            return record["returns"], record["ownership"], record["errors"]
        return (
            "The result declared for this name in include/seams.h.",
            "The caller owns any buffer it allocated. Library text is not freed by the caller.",
            STATUS if name.startswith("seams_") else "None.",
        )
    if lang == "cxx":
        record = cxx_records().get(name)
        if record:
            return record["returns"], record["ownership"], record["errors"]
        if name.isupper():
            return ("The constant value.", "No allocation.", "None.")
        return (
            "Zero on success. A non-zero value is the seams_status code.",
            "The caller owns every buffer it passed. The frame stays with its owner.",
            "A non-zero code. seams_last_error() holds the detail.",
        )
    return ("The value named in the purpose.", "The caller owns returned containers.", "A failing status becomes the language's error.")


def _field(block, label):
    match = re.search(r"^" + label + r": (.*)$", block, re.M)
    return match.group(1).strip() if match else ""


def _example(block, lang, name):
    current = _field(block, "Example").strip("`")
    if current:
        return current
    return EXAMPLES[lang]


HEADER_RENDERED = "The tour call is in `docs/examples/c_api_tour.c`."


def _complete_header_section(name, old, lang, record=None):
    """Add the missing contract lines to a section docs/capi_docs.py rendered from seams.h."""
    purpose = _field(old, "Purpose")
    parameters = _field(old, "Parameters")
    if record:
        returns, ownership, errors = record["returns"], record["ownership"], record["errors"]
    else:
        returns, ownership, errors = semantics(lang, name, purpose, parameters)
    body = old.strip("\n")
    for label, text in (("Ownership", ownership), ("Errors", errors)):
        if not re.search(r"^" + label + r": \S", body, re.M):
            body += "\n\n{}: {}".format(label, text[:1].upper() + text[1:])
    if not re.search(r"^Example: `", body, re.M):
        body += "\n\nExample: `{}`".format(_example(old, lang, name))
    return "### {}\n\n{}\n".format(name, body)


def render_section(name, old, lang, record):
    if lang == "c" and HEADER_RENDERED in old:
        return _complete_header_section(name, old, lang, record)
    purpose = _field(old, "Purpose")
    parameters = _field(old, "Parameters")
    units = _field(old, "Units")
    if purpose.startswith("Public name in") or not purpose:
        purpose = record["purpose"] if record else purpose or name
    if parameters.startswith("see the declaration") or not parameters:
        parameters = record["parameters"] if record else parameters or "none"
    if units.startswith("as documented on that declaration") or not units:
        units = record["units"] if record else units or "as named by the parameters"
    if record:
        returns, ownership, errors = record["returns"], record["ownership"], record["errors"]
    else:
        returns, ownership, errors = semantics(lang, name, purpose, parameters)
    if lang == "cxx":
        returns = _field(old, "Returns") or returns
        ownership = _field(old, "Ownership") or ownership
        errors = _field(old, "Errors") or errors
    example = _example(old, lang, name)
    if record and record.get("example") and not _field(old, "Example"):
        example = record["example"]
    known = ("Purpose", "Parameters", "Returns", "Ownership", "Errors", "Units", "Example")
    extra = []
    for line in old.splitlines():
        if line.startswith(tuple(item + ":" for item in known)) or not line.strip():
            continue
        extra.append(line)
    def cap(text):
        text = text.strip()
        if text and text[0].islower():
            text = text[0].upper() + text[1:]
        return text

    lines = [
        "### {}".format(name),
        "",
        "Purpose: {}".format(cap(purpose)),
        "",
        "Parameters: {}".format(parameters),
        "",
        "Returns: {}".format(cap(returns)),
        "",
        "Ownership: {}".format(cap(ownership)),
        "",
        "Errors: {}".format(cap(errors)),
        "",
        "Units: {}".format(units),
        "",
        "Example: `{}`".format(example),
        "",
    ]
    if extra:
        lines.extend(extra)
        lines.append("")
    return "\n".join(lines)


def sections_of(body):
    found = []
    chunks = re.split(r"^### ", body, flags=re.M)
    for chunk in chunks[1:]:
        line, _, rest = chunk.partition("\n")
        found.append((line.strip(), rest))
    return found


def records_for(lang):
    if lang == "c":
        return c_records()
    if lang == "cxx":
        return cxx_records()
    return {}


def public_names(lang):
    groups = {
        "c": check_c_symbols.c_header_symbols,
        "cxx": check_c_symbols.cxx_symbols,
    }
    return groups[lang]()


def rewritten_text(lang, text):
    start = text.find(BEGIN)
    stop = text.find(END)
    body = text[start:stop] if start >= 0 and stop > start else ""
    existing = sections_of(body)
    have = {}
    for name, rest in existing:
        have.setdefault(name, rest)
    records = records_for(lang)
    names = list(public_names(lang))
    if lang == "cxx":
        name_set = set(names)
        ordered = []
        seen = set()
        for name, _ in existing:
            if name in name_set and name not in seen:
                ordered.append(name)
                seen.add(name)
        for name in names:
            if name not in seen:
                ordered.append(name)
                seen.add(name)
    else:
        ordered = names
    blocks = ["", "## Symbol reference", ""]
    for name in ordered:
        record = records.get(name)
        blocks.append(render_section(name, have.get(name, ""), lang, record))
    wrapped = BEGIN + "\n" + "\n".join(blocks).rstrip() + "\n" + END
    if start >= 0 and stop > start:
        updated = text[:start] + wrapped + text[stop + len(END) :]
    else:
        if not text.endswith("\n"):
            text += "\n"
        updated = text + "\n" + wrapped + "\n"
    missing = [name for name in ordered if name not in have and name not in records and lang == "c"]
    return updated, len(ordered), missing


def rewrite(lang, api_dir=None, write=True):
    path = (Path(api_dir) if api_dir is not None else API) / PAGES[lang]
    original = path.read_text()
    updated, count, missing = rewritten_text(lang, original)
    if write:
        path.write_text(updated)
    return count, missing, original, updated


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument(
        "--api-dir",
        type=Path,
        default=API,
        help="Directory holding c.md and cpp.md (default: docs/api).",
    )
    parser.add_argument(
        "--check",
        action="store_true",
        help="Fail with a unified diff if rewriting would change c.md or cpp.md.",
    )
    args = parser.parse_args(argv)
    notes = []
    changed = False
    for lang in PAGES:
        count, missing, original, updated = rewrite(lang, api_dir=args.api_dir, write=not args.check)
        if args.check and updated != original:
            changed = True
            diff = difflib.unified_diff(
                original.splitlines(keepends=True),
                updated.splitlines(keepends=True),
                fromfile=str(args.api_dir / PAGES[lang]),
                tofile="generated/" + PAGES[lang],
            )
            sys.stderr.writelines(diff)
        notes.append("{} {}".format(lang, count))
        for name in missing:
            notes.append("  c fallback " + name)
    if changed:
        return 1
    print("\n".join(notes))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
