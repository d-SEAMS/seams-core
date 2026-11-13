#!/usr/bin/env python3
# Copyright (c) d-SEAMS developers
"""Write docs/api/c.md from include/seams.h and the tour's printed lines."""

import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
HEADER = ROOT / "include" / "seams.h"
TOUR = ROOT / "docs" / "examples" / "c_api_tour.c"
OUT = ROOT / "docs" / "api" / "c.md"

INTRO = """# C API

`seams.h` is ABI 1 of libseams. A handle is opaque. Coordinates are Angstrom.
A cell is nine numbers: `ox, oy, oz, lx, ly, lz, xy, xz, yz`. A zero tilt is orthorhombic.
`lx, ly, lz` are edge lengths.

`docs/examples/c_api_tour.c` calls every exported function on the in-tree fixtures.
`meson test docs_c_api_tour` compiles that file, runs it, and requires each printed
`capi` line below to match the run.

## Handle lifecycle

`seams_system_create` allocates an empty frame. The caller owns it and passes it to
`seams_system_destroy`. `NULL` is a valid argument to every destroy function.
`seams_system_set_atoms` stores points and a nine-number cell. Lengths are `size_t`.
`seams_system_set_frame` is the same store with an `int` count.
`seams_system_set_ortho` stores an orthorhombic cell. `seams_system_read` copies one
1-based frame from a trajectory. A negative type filter keeps every atom. A positive
filter keeps that species. `seams_system_count` is the number of points, or 0 for a
null handle. Neighbour, ring, cage, brine, and fingerprint handles are separate
objects. Destroying the frame does not free them.

One thread uses a handle. Another thread may use a different handle at the same time.
`seams_ql`, `seams_q4q6`, `seams_q6bar`, and `seams_steinhardt` share one
process-global scratch, so those four calls do not overlap.

## Status codes

`seams_status` is zero on success.

| code | value | `seams_status_string` |
| --- | --- | --- |
| `SEAMS_OK` | 0 | ok |
| `SEAMS_ERR_NULL` | 1 | null pointer |
| `SEAMS_ERR_RANGE` | 2 | value out of range |
| `SEAMS_ERR_BOX` | 3 | box is not a positive cell |
| `SEAMS_ERR_ARG` | 4 | bad argument |
| `SEAMS_ERR_STATE` | 5 | system has no frame |
| `SEAMS_ERR_NOMEM` | 6 | out of memory |
| `SEAMS_ERR_INTERNAL` | 7 | internal error |

## seams_last_error and seams_handle_error

`seams_last_error` is the detail of the last failing call on this thread. A later
successful call on the same thread clears it. The storage is thread-local.
`seams_handle_error` is the copy stored on the handle that failed. That copy stays
until the next call on that handle, or until destroy. Pass `NULL` and
`seams_handle_error` returns the thread-local string. The caller does not free either
pointer.

On an empty frame, `seams_chill_plus` returns `SEAMS_ERR_STATE`. The tour prints:

```
capi seams_last_error system has no frame
capi seams_handle_error no-frame
```

## Threads

`seams_set_num_threads(1)` selects one OpenMP thread. `seams_set_num_threads(0)`
restores the library default. A negative width returns `SEAMS_ERR_RANGE`.
`seams_get_num_threads` reports the width the next call will use.
`seams_set_num_threads` is process-global. Do not call it while another thread is
inside `seams_get_num_threads`.

## What the tour runs

The eight-atom cell is the diamond in `docs/examples/ice.c`. Every atom is cubic.
`input/traj/ice1h.lammpstrj` is the 128-oxygen ice Ih cell. `input/traj/genice_sI.lammpstrj`
is the 46-oxygen sI cage. `input/traj/exampleTraj.lammpstrj` is the mixed TIP4P dump
used for hydrogen bonds, sites, and domains. The cutoff is 3.5 angstroms unless a
line below says otherwise.

## Functions

"""


def symbols(text):
    found = set(re.findall(r"\bseams_[A-Za-z0-9_]+\b", text))
    found.update(re.findall(r"\bSEAMS_[A-Z0-9_]+\b", text))
    found.discard("SEAMS_API")
    return sorted(name for name in found if not name.endswith("_H_"))


def lines_for(name, printed):
    return [row for row in printed if row.split()[1] == name]


def main():
    header = HEADER.read_text()
    printed = [row for row in sys.stdin.read().splitlines() if row.startswith("capi ")]
    if not printed:
        print("tour output is empty", file=sys.stderr)
        return 1
    names = symbols(header)
    chunks = [INTRO]
    seen = []
    for row in printed:
        name = row.split()[1]
        if name in seen:
            chunks.append("```\n{}\n```\n".format(row))
            continue
        seen.append(name)
        chunks.append("### {}\n".format(name))
        chunks.append(
            "The call is in `docs/examples/c_api_tour.c`. "
            "The run prints:\n"
        )
        chunks.append("```\n{}\n```\n".format(row))
    chunks.append("<!-- API_SYMBOLS_BEGIN -->\n")
    chunks.append("## Symbol reference\n")
    for name in names:
        chunks.append("### {}\n".format(name))
        chunks.append("Purpose: Public name in include/seams.h.\n")
        chunks.append("Parameters: see the declaration in include/seams.h.\n")
        chunks.append("Units: as documented on that declaration.\n")
        chunks.append("Example: `docs/examples/c_api_tour.c`\n")
    chunks.append("<!-- API_SYMBOLS_END -->\n")
    chunks.append("The program those lines come from:\n")
    chunks.append("```{literalinclude} ../../../docs/examples/c_api_tour.c\n```\n")
    text = "\n".join(chunks)
    missing = [name for name in names if name not in text]
    if missing:
        print("generator missed " + " ".join(missing), file=sys.stderr)
        return 1
    OUT.write_text(text)
    print("wrote {} ({} capi lines)".format(OUT, len(printed)))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
