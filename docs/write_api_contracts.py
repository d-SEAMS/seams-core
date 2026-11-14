#!/usr/bin/env python3
# Copyright (c) d-SEAMS developers
"""Fill return, ownership, and error lines on the API reference pages.

Existing purpose, parameter, unit, and example lines stay. A placeholder
purpose is replaced with the header's brief. New public symbols gain a
section. Nothing outside the symbol block is removed.
"""

import re
from pathlib import Path

import check_c_symbols

ROOT = Path(__file__).resolve().parents[1]
API = ROOT / "docs" / "api"
BEGIN = "<!-- API_SYMBOLS_BEGIN -->"
END = "<!-- API_SYMBOLS_END -->"
PAGES = {
    "c": "c.md",
    "cxx": "cpp.md",
    "python": "python.md",
    "rust": "rust.md",
    "fortran": "fortran.md",
    "julia": "julia.md",
    "r": "r.md",
    "lua": "lua.md",
}
EXAMPLES = {
    "c": "docs/examples/c_api_tour.c",
    "cxx": "docs/examples/minimal_cpp.cpp",
    "python": "docs/examples/minimal_py.py",
    "rust": "docs/examples/minimal_rust.rs",
    "fortran": "docs/examples/minimal_fortran.f90",
    "julia": "docs/examples/minimal_julia.jl",
    "r": "docs/examples/minimal_r.R",
    "lua": "docs/examples/minimal_lua.lua",
}
FRAME_EXAMPLE = "docs/orgmode/reference/python-frame.org"
SEAMS_JL_EXAMPLE = "bindings/julia/Seams/test/runtests.jl"
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


def lua_records():
    text = (ROOT / "bindings" / "lua" / "dseams.lua").read_text()
    aliases = dict(re.findall(r"^dseams\.([A-Za-z0-9_]+)\s*=\s*dseams\.([A-Za-z0-9_]+)", text, re.M))
    hyphen = dict(re.findall(r'^dseams\["([^"]+)"\]\s*=\s*dseams\.([A-Za-z0-9_]+)', text, re.M))
    signatures = dict(re.findall(r"^function dseams\.([A-Za-z0-9_]+)\(([^)]*)\)", text, re.M))
    found = {}
    for name in check_c_symbols.lua_symbols():
        if name in aliases:
            purpose = "Alias of dseams.{}.".format(aliases[name])
            params = signatures.get(aliases[name], "none").replace(", ", "; ") or "none"
        elif name in hyphen:
            purpose = "Fennel spelling of dseams.{}.".format(hyphen[name])
            params = signatures.get(hyphen[name], "none").replace(", ", "; ") or "none"
        elif name in {"PUBLISHED", "OVITO", "CUTOFF", "MUTUAL", "SEEDED", "WATER_O"}:
            purpose = "Integer code published on the dseams table."
            params = "none"
        else:
            purpose = LUA_PURPOSE.get(name, "Lua binding for the matching C call.")
            params = signatures.get(name, "none").replace(", ", "; ") or "none"
        if name in {"PUBLISHED", "OVITO", "CUTOFF", "MUTUAL", "SEEDED", "WATER_O"}:
            returns, ownership, errors = (
                "The integer code.",
                "No allocation.",
                "None.",
            )
        elif name == "cloud.close":
            returns, ownership, errors = (
                "Nothing.",
                "Frees the seams_system inside the cloud. A second call is a no-op.",
                "None.",
            )
        elif name == "core":
            returns, ownership, errors = (
                "The dseams_core module table.",
                "Lua owns the module. The caller does not free C pointers from it by hand; the wrappers do.",
                "A failing C call raised by a wrapper is a Lua error whose text is seams_last_error.",
            )
        else:
            returns = "The Lua table or userdata described for this call."
            ownership = (
                "A system userdata is freed by cloud:close(). "
                "Lua 5.4 also frees it at the end of a <close> scope. "
                "Returned tables are owned by Lua."
            )
            errors = "A status other than SEAMS_OK is a Lua error whose text is seams_last_error."
        found[name] = {
            "purpose": _sentence(purpose),
            "parameters": params if params else "none",
            "returns": returns,
            "ownership": ownership,
            "errors": errors,
            "units": "dimensionless code" if name.isupper() else "cutoff and coordinates in angstrom when the call takes them",
            "example": EXAMPLES["lua"],
        }
    return found


def python_extra():
    frame = {
        "Frame.set_frame": (
            "Store positions and a cell on a dseams.Frame.",
            "positions; box",
            "Nothing.",
            "The frame copies the coordinates into the seams_system it owns.",
            "ValueError when box is not 3 lengths or 9 LAMMPS numbers. RuntimeError with seams_last_error text when the C call fails.",
        ),
        "Frame.count": (
            "Return the number of points stored on a dseams.Frame.",
            "self",
            "The atom count.",
            "No allocation.",
            "None.",
        ),
        "Frame.chill_plus": (
            "CHILL+ labels and the nine population counts for the atoms on a dseams.Frame.",
            "self; cutoff; graph; rule",
            "A pair of NumPy arrays, labels and counts.",
            "The caller owns the returned arrays. The frame still owns the seams_system.",
            "RuntimeError whose text is seams_last_error.",
        ),
        "Frame.q4q6": (
            "Local q4 and q6 on the cutoff shell of a dseams.Frame.",
            "self; cutoff",
            "A pair of NumPy arrays, q4 and q6.",
            "The caller owns the returned arrays.",
            "RuntimeError whose text is seams_last_error.",
        ),
        "Frame.neighbours": (
            "CSR neighbour list of a dseams.Frame, self entries omitted.",
            "self; cutoff; graph",
            "A pair of NumPy arrays, offsets and columns.",
            "The caller owns the returned arrays. The temporary C neighbour list is freed inside the call.",
            "RuntimeError whose text is seams_last_error.",
        ),
        "Frame.rings": (
            "Ring sizes, symmetry labels, a size histogram, and the six-ring fraction.",
            "self; cutoff; max_depth",
            "Sizes, symmetry, histogram, and the six-ring fraction.",
            "The caller owns the returned arrays. The C ring set is freed inside the call.",
            "RuntimeError whose text is seams_last_error.",
        ),
    }
    user = {
        "userframe.IceCounts": (
            "CHILL+ population counts keyed by label name.",
            "mapping of name to count",
            "An IceCounts dictionary.",
            "The caller owns the dictionary.",
            "None.",
        ),
        "userframe.CageScore": (
            "Boolean membership of the hexagonal and double-diamond cages.",
            "hc; ddc",
            "A CageScore.",
            "The caller owns the object.",
            "None.",
        ),
        "userframe.CageScore.n_ih": (
            "Number of atoms in the hexagonal cage.",
            "self",
            "The count.",
            "No allocation.",
            "None.",
        ),
        "userframe.CageScore.n_ic": (
            "Number of atoms in the double-diamond cage.",
            "self",
            "The count.",
            "No allocation.",
            "None.",
        ),
        "userframe.CageScore.n_water": (
            "Number of atoms in neither cage.",
            "self",
            "The count.",
            "No allocation.",
            "None.",
        ),
        "userframe.CageScore.n_both": (
            "Number of atoms in both cages.",
            "self",
            "The count.",
            "No allocation.",
            "None.",
        ),
        "userframe.Frame": (
            "One configuration whose classifiers call libseams_lang.",
            "ptr; symbols; cutoff",
            "A Frame.",
            "The object owns the seams_frame pointer. close() or deletion calls seams_frame_free.",
            "RuntimeError when ptr is null.",
        ),
        "userframe.Frame.close": (
            "Free the libseams_lang frame.",
            "self",
            "Nothing.",
            "Frees the pointer this Frame owns. A second call is a no-op.",
            "None.",
        ),
        "userframe.Frame.n_atoms": (
            "Number of atoms on the frame.",
            "self",
            "The count from seams_frame_nop.",
            "No allocation.",
            "None.",
        ),
        "userframe.Frame.positions": (
            "Copy of the coordinates stored on the frame.",
            "self",
            "An (n, 3) NumPy array in angstrom.",
            "The caller owns the copy. The frame keeps its own buffer.",
            "None.",
        ),
        "userframe.Frame.box": (
            "Cell stored on the frame.",
            "self",
            "Three edge lengths, or nine LAMMPS numbers when a tilt is set.",
            "The caller owns the returned array.",
            "None.",
        ),
        "userframe.Frame.chill_plus": (
            "Published CHILL+ histogram for this frame.",
            "self; graph",
            "An IceCounts mapping of label name to population.",
            "The caller owns the dictionary. Labels are also stored on the frame.",
            "RuntimeError whose text is seams_last_error.",
        ),
        "userframe.Frame.cages": (
            "Hexagonal and double-diamond membership from seams_frame_cages.",
            "self; seeded; k; ring_adjacent",
            "A CageScore.",
            "The caller owns the score. The frame keeps the same object.",
            "RuntimeError whose text is seams_last_error.",
        ),
        "userframe.Frame.cages_by_signature": (
            "Closed cages of one ring census.",
            "self; signature; max_ring_size",
            "A list whose length is the cage count. n_atoms is the vertex count.",
            "The caller owns the list. The C call does not return one vertex list per cage.",
            "RuntimeError whose text is seams_last_error.",
        ),
        "userframe.Frame.steinhardt": (
            "Local q_l and the Lechner-Dellago average from seams_frame_steinhardt.",
            "self; order_l",
            "A dictionary with ql and ql_bar lists.",
            "The caller owns the dictionary.",
            "RuntimeError whose text is seams_last_error.",
        ),
        "userframe.Frame.ring_census": (
            "Primitive rings by size from seams_frame_rings.",
            "self; max_size",
            "A dictionary of size to count, empty sizes omitted.",
            "The caller owns the dictionary.",
            "RuntimeError whose text is seams_last_error.",
        ),
        "userframe.Frame.rdf_2d": (
            "In-plane g(r) from seams_frame_rdf2d.",
            "self; output_dir; cutoff; binwidth",
            "Radius and g(r) lists.",
            "The caller owns the lists. The call also writes under output_dir.",
            "RuntimeError whose text is seams_last_error.",
        ),
        "userframe.Frame.monolayer_rings": (
            "Quasi-2D ring counts from seams_frame_monolayer.",
            "self; output_dir; sheet_area; max_depth",
            "A dictionary of size to count and coverage.",
            "The caller owns the dictionary. The call writes topoMonolayer under output_dir.",
            "RuntimeError whose text is seams_last_error.",
        ),
        "userframe.Frame.find_prisms": (
            "Prism census from seams_frame_prisms.",
            "self; output_dir; max_depth; shape_matching",
            "A list of size, perfect, and imperfect rows.",
            "The caller owns the list.",
            "RuntimeError whose text is seams_last_error.",
        ),
        "userframe.Frame.to_ase": (
            "ASE Atoms for this frame.",
            "self",
            "An ASE Atoms object. Ice labels and cage flags are arrays when those calls have run.",
            "The caller owns the Atoms object.",
            "ImportError when ASE is not installed.",
        ),
        "userframe.Frame.to_solvis": (
            "solvis System built from to_ase.",
            "self; expand_box",
            "A solvis System.",
            "The caller owns the System.",
            "ImportError when solvis-tools is not installed.",
        ),
        "userframe.read": (
            "Read one frame with seams_frame_read. pydseams.read is this function.",
            "filename; frame; atom_type; all_atoms",
            "A userframe.Frame.",
            "The caller owns the Frame and frees it with close() or by dropping it.",
            "RuntimeError whose text is seams_last_error.",
        ),
        "userframe.from_ase": (
            "Frame for one ASE species. Coordinates go through seams_frame_from_arrays.",
            "atoms; select; cutoff; bonded",
            "A userframe.Frame.",
            "The caller owns the Frame.",
            "TypeError when atoms has no positions. ValueError when the selection matches nothing.",
        ),
        "userframe.to_solvis": (
            "Call to_solvis on a frame.",
            "frame; expand_box",
            "The solvis System returned by the frame.",
            "The caller owns the System.",
            "ImportError when solvis-tools is not installed.",
        ),
    }
    rows = {}
    for name, fields in {**frame, **user}.items():
        rows[name] = {
            "purpose": _sentence(fields[0]),
            "parameters": fields[1],
            "returns": _sentence(fields[2]),
            "ownership": _sentence(fields[3]),
            "errors": _sentence(fields[4]),
            "units": "coordinates and cutoff in angstrom",
            "example": FRAME_EXAMPLE if name.startswith("userframe.") else EXAMPLES["python"],
        }
    return rows


def julia_extra():
    rows = {
        "ice_ih_path": (
            "Path of the committed Ice Ih dump, input/traj/ice1h.lammpstrj.",
            "none",
            "The absolute path string.",
            "Julia owns the string.",
            "None.",
            EXAMPLES["julia"],
        ),
        "Seams.System": (
            "Empty system from seams_system_create. Positions are added with set_ortho!.",
            "none",
            "A System.",
            "The finalizer calls seams_system_destroy.",
            "error when seams_system_create returns NULL.",
            SEAMS_JL_EXAMPLE,
        ),
        "Seams.version": (
            "Engine version string from seams_version_string.",
            "none",
            "A String.",
            "Julia owns the string. The C pointer is not freed.",
            "None.",
            SEAMS_JL_EXAMPLE,
        ),
        "Seams.set_num_threads": (
            "Request n worker threads.",
            "n",
            "Nothing.",
            "No allocation.",
            "error with seams_last_error text when the status is not SEAMS_OK.",
            SEAMS_JL_EXAMPLE,
        ),
        "Seams.num_threads": (
            "Worker threads the engine is using.",
            "none",
            "An Int.",
            "No allocation.",
            "None.",
            SEAMS_JL_EXAMPLE,
        ),
        "Seams.set_ortho!": (
            "Store an orthogonal cell. xyz is 3 by n.",
            "sys; xyz; lengths",
            "The same System.",
            "The System still owns the C frame.",
            "error when xyz is not 3 by n, or when the C status is not SEAMS_OK.",
            SEAMS_JL_EXAMPLE,
        ),
        "Seams.chill_plus": (
            "CHILL+ labels and the nine-bin census.",
            "sys; cutoff; graph; rule",
            "A pair of vectors, labels and counts.",
            "The caller owns the vectors. The System still owns the frame.",
            "error with seams_last_error text when the status is not SEAMS_OK.",
            SEAMS_JL_EXAMPLE,
        ),
        "Seams.q4q6": (
            "Local q4 and q6.",
            "sys; cutoff",
            "A pair of vectors, q4 and q6.",
            "The caller owns the vectors.",
            "error with seams_last_error text when the status is not SEAMS_OK.",
            SEAMS_JL_EXAMPLE,
        ),
        "Seams.six_rings": (
            "Number of primitive six-rings and their fraction.",
            "sys; cutoff",
            "A pair, the six-ring count and the fraction.",
            "The temporary ring set is freed with seams_rings_destroy before the function returns.",
            "error with seams_last_error text when the status is not SEAMS_OK.",
            SEAMS_JL_EXAMPLE,
        ),
        "Seams.PUBLISHED": (
            "CHILL+ rule code of the 2015 paper.",
            "none",
            "The integer 0.",
            "No allocation.",
            "None.",
            SEAMS_JL_EXAMPLE,
        ),
        "Seams.OVITO": (
            "CHILL+ rule code used by OVITO.",
            "none",
            "The integer 1.",
            "No allocation.",
            "None.",
            SEAMS_JL_EXAMPLE,
        ),
        "Seams.CUTOFF": (
            "Graph code of the cutoff neighbour list.",
            "none",
            "The integer 0.",
            "No allocation.",
            "None.",
            SEAMS_JL_EXAMPLE,
        ),
        "Seams.MUTUAL": (
            "Graph code of the mutual four-nearest list.",
            "none",
            "The integer 1.",
            "No allocation.",
            "None.",
            SEAMS_JL_EXAMPLE,
        ),
    }
    out = {}
    for name, fields in rows.items():
        out[name] = {
            "purpose": _sentence(fields[0]),
            "parameters": fields[1],
            "returns": _sentence(fields[2]),
            "ownership": _sentence(fields[3]),
            "errors": _sentence(fields[4]),
            "units": "cutoff and coordinates in angstrom when the call takes them",
            "example": fields[5],
        }
    return out


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
        if name.endswith("_free") or name.endswith("_destroy"):
            return (
                "Nothing.",
                "Frees the caller-owned frame. A null pointer is a no-op.",
                "None.",
            )
        if name in {"seams_frame_read", "seams_frame_from_arrays", "seams_frame_from_species"}:
            return (
                "A new seams_frame, or NULL when the call cannot build one.",
                "The caller owns the frame and frees it with seams_frame_free.",
                "NULL. seams_last_error() names the reason.",
            )
        if name in {"seams_frame_positions", "seams_frame_types", "seams_frame_ids"}:
            return (
                "A pointer into the frame.",
                "The frame owns the memory. It stays valid until seams_frame_free or a call that replaces the coordinates.",
                "NULL when the frame is null.",
            )
        if name in {"seams_last_error", "seams_lang_version"}:
            return (
                "A pointer to library text.",
                "The library owns the string. The caller does not free it.",
                "None.",
            )
        if name in {"seams_frame", "seams_phase_times", "seams_prism_row"}:
            return ("The type itself.", "No allocation. This is a type.", "None.")
        return (
            "Zero on success. A non-zero value is the seams_status code.",
            "The caller owns every buffer it passed. The frame stays with its owner.",
            "A non-zero code. seams_last_error() holds the detail.",
        )
    if lang == "python":
        extra = python_extra().get(name)
        if extra:
            return extra["returns"], extra["ownership"], extra["errors"]
        if name in {"CUTOFF", "MUTUAL", "OVITO", "PUBLISHED"}:
            return ("The integer code.", "No allocation.", "None.")
        if name == "Frame":
            return (
                "A Frame object.",
                "The object owns its seams_system and frees it when the object is destroyed.",
                "None on construction. Later methods raise RuntimeError with seams_last_error text.",
            )
        return (
            "The arrays or counts named in the purpose. NumPy arrays are owned by the caller.",
            "A Frame owns its seams_system and frees it on destruction. Input arrays are read in place when they are already C-contiguous.",
            "RuntimeError with seams_last_error text. ValueError when a box is not 3 lengths or 9 LAMMPS numbers.",
        )
    if lang == "rust":
        if name in {"System", "Neighbours", "Rings", "Cages", "Fingerprint"}:
            return (
                "The owned handle.",
                "Drop frees the C object (seams_system_destroy or the matching destroy function).",
                "Constructors return Err(Error) when the C call fails. Error.code is the seams_status.",
            )
        if name == "Error":
            return (
                "The status code and the detail string.",
                "The caller owns the String.",
                "This is the error type. code is the seams_status value and message is seams_handle_error or seams_last_error.",
            )
        if name in {"version", "count", "num_threads", "nnz"}:
            return (
                "The value named in the purpose.",
                "A returned String is owned by the caller. No C handle is allocated.",
                "None for the integer queries. version copies library text.",
            )
        return (
            "Result of the value, or Err(Error).",
            "Returned Vec and String values are owned by the caller. Handles returned as structs are freed by Drop.",
            "Error.code is the seams_status. Error.message is seams_handle_error when the call had a handle, otherwise seams_last_error.",
        )
    if lang == "fortran":
        if name == "c_chars":
            return (
                "An allocatable C string, including the trailing null.",
                "The caller owns the allocation.",
                "An allocation failure stops the program.",
            )
        if name.startswith("SEAMS_"):
            return ("The integer parameter.", "No allocation.", "None.")
        if name.endswith("_destroy"):
            return (
                "Nothing.",
                "Frees the c_ptr the caller owns.",
                "None.",
            )
        return (
            "integer(c_int) status. SEAMS_OK is zero.",
            "A c_ptr from seams_system_create is freed with seams_system_destroy. Allocatable character buffers from c_chars are owned by the caller.",
            "A non-zero seams_status. seams_last_error and seams_handle_error hold the detail.",
        )
    if lang == "julia":
        extra = julia_extra().get(name)
        if extra:
            return extra["returns"], extra["ownership"], extra["errors"]
        if name == "Frame":
            return (
                "A Frame.",
                "The finalizer calls seams_frame_free.",
                "error when the pointer is null.",
            )
        return (
            "The Julia value named in the purpose.",
            "Returned arrays are owned by Julia. A Frame is freed by its finalizer.",
            "error() with the text of seams_last_error when the status is not zero.",
        )
    if lang == "r":
        return (
            "The data frame, vector, or count named in the purpose.",
            "R owns the returned vectors. A dseams_frame finalizer frees the C frame.",
            "The binding stops with the text of seams_last_error when the status is not SEAMS_OK.",
        )
    if lang == "lua":
        record = lua_records().get(name)
        if record:
            return record["returns"], record["ownership"], record["errors"]
    return ("The value named in the purpose.", "The caller owns returned containers.", "A failing status becomes the language's error.")


def _field(block, label):
    match = re.search(r"^" + label + r": (.*)$", block, re.M)
    return match.group(1).strip() if match else ""


def _example(block, lang, name):
    current = _field(block, "Example").strip("`")
    if current:
        return current
    if lang == "python" and name.startswith("userframe."):
        return FRAME_EXAMPLE
    if lang == "julia" and (name.startswith("Seams.") or name == "ice_ih_path"):
        record = julia_extra().get(name)
        if record:
            return record["example"]
    return EXAMPLES[lang]


def render_section(name, old, lang, record):
    purpose = _field(old, "Purpose")
    parameters = _field(old, "Parameters")
    units = _field(old, "Units")
    if purpose.startswith("Public name in") or not purpose:
        purpose = record["purpose"] if record else purpose or name
    if parameters.startswith("see the declaration") or not parameters:
        parameters = record["parameters"] if record else parameters or "none"
    if units.startswith("as documented on that declaration") or not units:
        units = record["units"] if record else units or "as named by the parameters"
    returns, ownership, errors = semantics(lang, name, purpose, parameters)
    if record and lang in {"lua", "python", "julia"} and name in (record if isinstance(record, dict) else {}):
        pass
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
    if lang == "lua":
        return lua_records()
    if lang == "python":
        return python_extra()
    if lang == "julia":
        return julia_extra()
    return {}


def public_names(lang):
    groups = {
        "c": check_c_symbols.c_header_symbols,
        "cxx": check_c_symbols.cxx_symbols,
        "python": check_c_symbols.python_symbols,
        "rust": check_c_symbols.rust_symbols,
        "fortran": check_c_symbols.fortran_symbols,
        "julia": check_c_symbols.julia_symbols,
        "r": check_c_symbols.r_symbols,
        "lua": check_c_symbols.lua_symbols,
    }
    return groups[lang]()


def rewrite(lang):
    path = API / PAGES[lang]
    text = path.read_text()
    start = text.find(BEGIN)
    stop = text.find(END)
    body = text[start:stop] if start >= 0 and stop > start else ""
    have = {name: rest for name, rest in sections_of(body)}
    records = records_for(lang)
    names = list(public_names(lang))
    blocks = ["", "## Symbol reference", ""]
    for name in names:
        record = records.get(name)
        blocks.append(render_section(name, have.get(name, ""), lang, record))
    wrapped = BEGIN + "\n" + "\n".join(blocks).rstrip() + "\n" + END + "\n"
    if start >= 0 and stop > start:
        text = text[:start] + wrapped + text[stop + len(END) :]
        if not text[stop + len(END) :].startswith("\n"):
            # END was consumed with the following newline in the slice above when we cut at stop.
            pass
        # The slice text[stop + len(END):] drops a single newline that followed END if we want to keep it.
        # wrapped already ends with newline. If the old END was followed by newline, we now have a double
        # newline only when the original had END\\n and we didn't include it. That's fine.
    else:
        if not text.endswith("\n"):
            text += "\n"
        text += "\n" + wrapped
    # Collapse the accidental gap from an END that already included a trailing newline in the source.
    text = text.replace(END + "\n\n\n", END + "\n\n")
    path.write_text(text)
    missing = [name for name in names if name not in have and name not in records and lang == "c"]
    return len(names), missing


def main():
    notes = []
    for lang in PAGES:
        count, missing = rewrite(lang)
        notes.append("{} {}".format(lang, count))
        for name in missing:
            notes.append("  c fallback " + name)
    print("\n".join(notes))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
