<!-- https://docs.dseams.info/explanation/architecture.md -->

# Architecture {#architecture}

d-SEAMS is one C++ engine, `libyodaLib`. Callers link
`libseams`, `libseams_lang`, or both. `libseams` (`include/seams.h`) is
the C API. `libseams_lang` (`bindings/cxx/seams_lang.h`) is the C ABI.
Each library links `libyodaLib`. `libseams_lang` is a separate library
from `libseams`. The `seams` command is the command-line front on the
same engine.

The engine reads the trajectory. XYZ and LAMMPS dumps are built in.
`.pdb`, `.gro`, and `.dcd` need chemfiles, and `.con` needs readcon.
Without chemfiles the command refuses `.pdb`, `.gro`, and `.dcd`.
Without readcon it refuses `.con`. Every other suffix is read as a
LAMMPS dump.

Optional drivers for LAMMPS, GROMACS, OpenMM, and PLUMED live in this
repository.

Building the engine and the `seams` command takes a C and C++
toolchain. Python, Lua, R, Julia, Rust, and Fortran are not required
for that build.

## Language packages {#language-packages}

Each language is its own package. The package calls this engine. It
does not carry a second copy of a classifier. A package release carries
the engine version it was built against, so the `pydseamslib` package
calls that engine.

| language | repository | package | calls |
| --- | --- | --- | --- |
| Python | PydSEAMSlib | pydseamslib | libseams (import dseams) and libseams_lang (import pydseams) |
| Lua, Fennel | yodaStruct | Git source, require("dseams") | libseams |
| R | dseams-r | dseams | libseams_lang |
| Julia | Dseams.jl | Dseams, Seams | Dseams calls libseams_lang. Seams calls libseams. |
| Rust | dseams-rs | dseams-sys, dseams | libseams |
| Fortran | dseams-fortran | fpm Git dependency dseams | libseams |

[Use a language binding](../user/howto/bindings.md) has the install
commands.

## Why the name is yoda {#why-the-name-is-yoda}

The 1.x executable was `yodaStruct`. This repository builds `seams`.
It does not build `yodaStruct`. The compiled Python module is
`pydseams.yoda`. `_core` and `cyoda` are aliases of `yoda`. Lua loads
`require("dseams")`. `require("yoda")` is that table.

## Cutoff list {#cutoff-list}

The cutoff list is every pair inside the distance cutoff.
[linkcell](https://github.com/d-SEAMS/linkcell) builds it, and there
is no other builder. `lc_pairs_within_csr` returns the rows as CSR,
built on linkcell’s threads and the same bytes at every width;
seams-core maps them to atom IDs or cloud indices. That covers the
full and half lists, orthorhombic and restricted triclinic boxes, and
cutoffs past half the box. This repository asks for linkcell 0.6.0 or
newer, and `subprojects/linkcell.wrap` pins `d-SEAMS/linkcell` main.
linkcell’s release number is its own, and a seams-core release tag does
not name it. `libyodaLib` links one copy and exports the `lc_` symbols.

A cutoff that is not finite and positive throws
`std::invalid_argument`, and a coordinate that is not finite throws
`linkcell::Error`. The C API returns nonzero with the message in
`seams_c_error()`. `seams --features` prints the compile-time backends;
`linkcell pairs_within` is always `enabled`.

The bond graph is a separate choice. For `seams cages`, `--graph` is
`cutoff`, `knn`, `knn-union`, or `seeded`, and `seeded` is the default
for that command. The cutoff-list builder does not select the graph.
Which graph the six-rings are found on, and why that is a different
question from CHILL+, is [CHILL+ versus cages](chill-plus-vs-cages.md).

## Device paths {#device-paths}

Two optional device builds do different work.

OpenMP target offload is the Meson option `with_openmp_offload`. The
default is disabled. The build covers Steinhardt $Y_{lm}$ and the
cage score, as hop-bound six-rings and cage affiliation, on the
neighbour list the host built. `SEAMS_OFFLOAD` requests that path.
The default is off, and `0` keeps the work on the host. `seams cages`
with `--graph seeded`, the default, stays on the host either way.
`seams chill` and `seams chill-plus` run on the host.

gpulite is the Meson option `with_gpulite`. When that library is
linked, the engine provides `gpu::analyzeResident`, a CUDA batch of
frames that stays on the device. The batch is the cage score: linkcell
device k-nearest, the mutual four-nearest graph, primitive six-rings,
and HC/DDC labels. It needs that device k-nearest. Without it, the
call reports that the path is not built. The batch returns cage
labels. It does not compute CHILL, CHILL+, or $q_{lm}$. The `seams`
subcommands do not call `gpu::analyzeResident`.

## Runtime settings {#runtime-settings}

Defaults are compiled into the binary. A dotenv file (`SEAMS_CONFIG`,
or `./seams.env` when that file is present) fills keys that are still
unset. The process environment overrides the file. A flag on the
`seams` command overrides the value for the setting that flag names.
`seams --print-config` prints the resolved values.

The keys choose frames, workers, the default graph, and whether
offload was requested. The subcommand (`chill`, `cages`, `rings`, and
the others) is not one of the keys. Compile-time Meson options are not
in the file. [The command reference](../reference/cli.md) lists
the keys and the defaults.

## Two classifiers {#two-classifiers}

`libyodaLib` holds both ice classifiers. CHILL and CHILL+ label a
molecule from its four nearest neighbours. The cage score labels
hexagonal-cage and double-diamond-cage membership. The two results do
not imply each other. The argument, and the counts on the mixed TIP4P
dump, are on [CHILL+ versus cages](chill-plus-vs-cages.md).
