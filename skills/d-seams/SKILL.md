---
name: d-seams
description: "Classify ice and hydrate cages in a molecular-dynamics trajectory with d-SEAMS. Use when installing the engine, running the seams CLI, calling libseams (seams.h) or pydseamslib (import dseams), or debugging CHILL+, cages, primitive rings, or Steinhardt q_l. Keywords: molecular dynamics, water, ice Ih, ice Ic, clathrate, cutoff, ABI."
license: MIT
compatibility: "The engine builds with Meson. Python is the package pydseamslib. Lua, R, Julia, Rust, and Fortran bindings live in their own repositories."
metadata:
  author: d-SEAMS
  version: "2.14.3"
---

# Use d-SEAMS

d-SEAMS reads a trajectory and names the ice. This repository is the
C++ engine (`libyodaLib`), the C API `libseams`, the C ABI `libseams_lang`,
and the `seams` command. Cite Goswami, Goswami, and Singh, J. Chem. Inf.
Model. 2020, https://doi.org/10.1021/acs.jcim.0c00031.

The manual is https://docs.dseams.info. The page index is
https://docs.dseams.info/llms.txt. Each HTML page has a Markdown twin at
the same path with `.md`.

## Install

From a source checkout, with Eigen 3.4 available:

```bash
meson setup build && meson compile -C build
./build/seams --version
```

`pixi run setup && pixi run build` uses the locked toolchain and writes
`./bbdir/seams`. `nix build` writes `./result/bin/seams`.

Python is a separate package. Building this repository does not install it.

```bash
pip install pydseamslib
```

`import dseams` is the analysis API. `import pydseams` is the lower-level
module from the same package.

## CLI

```bash
seams read FILE
seams chill-plus FILE --cutoff 3.5 --type 2
seams cages FILE --type 2
seams rings FILE --type 2
seams steinhardt FILE --type 2 --l 6
seams --version
seams --features
seams --help
```

The default neighbour cutoff is 3.5 angstrom. `--type` selects the species.
TIP4P oxygen in the bundled mixed dump is type 2. The mW and cubic-ice
dumps are type 1. `nop 0` means the selected type matched no atoms.

`seams chill-plus` uses the cutoff neighbour list. `seams cages` defaults
to `--graph seeded` (mutual four-nearest seeds, then union completions).
Cage membership and the CHILL+ star answer different questions. A mixed
frame can carry interfacial labels and still contain no complete HC or DDC.

Further commands include `chill`, `rdf`, `cn`, `hbonds`, `pairs`,
`density-z`, `domains`, `prisms`, `fingerprint`, `ions`, and `symmetry`.

Runtime knobs (`SEAMS_CUTOFF`, `SEAMS_TYPE`, `SEAMS_JOBS`, `SEAMS_GRAPH`,
and the rest of that table) resolve in this order: the built-in default,
then `SEAMS_CONFIG` or `./seams.env`, then the process environment, then
the flag. `seams --print-config` prints the resolved table.

## Python

```python
import dseams
a = 2.75 * 4 / 3**0.5
sites = ((0, 0, 0), (0.5, 0.5, 0), (0.5, 0, 0.5), (0, 0.5, 0.5),
         (0.25, 0.25, 0.25), (0.75, 0.75, 0.25), (0.75, 0.25, 0.75),
         (0.25, 0.75, 0.75))
xyz = [[x * a, y * a, z * a] for x, y, z in sites]
labels, counts = dseams.chill_plus(xyz, [a, a, a])
assert counts[0] == 8  # cubic
```

Coordinates are Cartesian, in angstrom. Pass orthorhombic lengths, or the
LAMMPS nine-value cell when the call takes one.

## C API

`include/seams.h` is ABI 1. Handles are opaque. A call that can fail
returns a `seams_status`. Zero is success. The detail string is
`seams_last_error_message()`, and that string is thread-local. The library
does not abort the process, and a C++ exception does not cross this
boundary.

```c
#include <seams.h>

seams_system *sys = seams_system_create();
seams_system_set_ortho(sys, xyz, n, 0, 0, 0, lx, ly, lz, NULL, NULL);
seams_chill_plus(sys, 3.5, 0, SEAMS_GRAPH_CUTOFF, SEAMS_CHILL_PUBLISHED,
                 labels, counts);
seams_system_destroy(sys);
```

Destroy every handle you create. `destroy` accepts NULL. A second destroy
of the same pointer records "double destroy" and frees nothing. Outputs
are caller buffers. A NULL output buffer with a non-NULL length pointer
queries the required count and returns success. A short buffer returns
`SEAMS_ERR_RANGE` and writes nothing into that buffer. Lengths are
`size_t`.

Distinct handles may be used from distinct threads at the same time. Two
non-const calls on one handle must not overlap. `seams_set_num_threads`
is process-global.

## Where the other languages live

| Language | Repository | Package |
| --- | --- | --- |
| Python | d-SEAMS/PydSEAMSlib | `pydseamslib` |
| Lua, Fennel | d-SEAMS/yodaStruct | Git source, `require("dseams")` |
| R | d-SEAMS/dseams-r | `dseams` |
| Julia | d-SEAMS/Dseams.jl | `Dseams` |
| Rust | d-SEAMS/dseams-rs | `dseams-sys`, `dseams` |
| Fortran | d-SEAMS/dseams-fortran | fpm `dseams` |

Land an engine change here first. The binding then moves its pin. Binding
source does not belong in this repository. The contributor instructions are the contract
for editing the tree.

## Pitfalls

- The 1.x command was `yodaStruct -c lua_inputs/config.yml`. The command in this tree is `seams`.
- A sheared LAMMPS dump already carries its tilt factors. Keep that cell. Replacing it with the diagonal of the bound spans drops pairs that cross a tilted face.
- conda-forge chemfiles 0.10.4 has 3 of the 37 VMD molfile readers and fails an mmCIF read of an atom with no residue. The Meson build probes it with `cmake/chemfiles_mmcif_probe.cpp` and vendors a patched chemfiles when the probe fails.
- `.con` needs readcon-core. The Nix flake does not fetch that wrap.
- `--jobs` needs an OpenMP build. `seams --features` reports whether OpenMP is on.
- `require("dseams")` needs the module built from a yodaStruct checkout.
- A bug report needs `seams --version`, `seams --features`, the command line, and a small dump.
