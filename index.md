<!-- https://docs.dseams.info/index.md -->

# `d-SEAMS` {#d-seams}

**Author**

[d-SEAMS core team](https://dseams.info)

Load a LAMMPS frame and print nop, frame, and box.

CHILL+ four-neighbour star: cubic, hexagonal, interfacial, clathrate, interClathrate, or water.

Cage membership: HC is ice Ih, DDC is ice Ic, else water.

## Overview {#overview}

d-SEAMS 2.14 reads a molecular-dynamics trajectory and names the ice in it.
The name is Deferred Structural Elucidation Analysis for Molecular
Simulations. Release 2.14.0 of this repository is the C++ engine
`libyodaLib`, the C API `libseams` (`include/seams.h`), the header-only
C++ API `include/seams.hpp` (it links `libseams`), the C ABI
`libseams_lang` (`bindings/cxx/seams_lang.h`, used by the Python, R, and
Julia fronts), the `seams` command, and the LAMMPS, GROMACS, OpenMM,
PLUMED, and in-situ drivers. `gmx-seams` reads the trajectory `gmx mdrun` wrote.
The standalone PLUMED action (`dseams-plumed` v0.4.1) is indexed in
[PLUMED-NEST egg 26.018](https://www.plumed-nest.org/eggs/26/018/).
`seams --version` prints `2.14.0`.

See the [PLUMED](engines/plumed.md) page for the plugin.

CHILL+ labels an oxygen from its four-neighbour star: cubic, hexagonal,
interfacial, clathrate, interClathrate, or water. Cage membership is a
separate score. A molecule in a hexagonal cage (HC) is ice Ih. A molecule
in a double-diamond cage (DDC) is ice Ic. A molecule in neither is water.
Primitive rings and Steinhardt $q_l$ are computed in the same engine.
[CHILL+ versus cages](explanation/chill-plus-vs-cages.md) records
where the two scores disagree.

A source build from this tree, when the compilers are already on `PATH`,
is `meson setup build && meson compile -C build`.
pixi is `pixi run setup && pixi run build`.
[Install](getting-started/install.md) has each method.

The ice Ih frame shipped in this tree is `input/traj/ice1h.lammpstrj`,
128 molecules of LAMMPS type 1.

```bash
seams read input/traj/ice1h.lammpstrj
seams chill-plus input/traj/ice1h.lammpstrj --cutoff 3.5
```

`seams read` prints `nop 128 frame 1 box 15.6457 14.7071 18.073`.
`seams chill-plus` prints `nop 128 hexagonal 128`.
[The quick start](getting-started/quickstart.md) classifies an eight-atom cubic
cell from C and from Python.
[From a trajectory to the ice front](user/howto/from-trajectory.md)
adds the ring census on this frame and the ice front on
`input/traj/interface_f0.lammpstrj`.

Each other language is its own repository.
[Use a language binding](user/howto/bindings.md) names the package
and the library each one links.

- Python: `pydseamslib` (`import dseams`, `import pydseams`),
  with `pip install 'pydseamslib==2.14.0'` from
  [PydSEAMSlib](https://github.com/d-SEAMS/PydSEAMSlib)
- Lua and Fennel: `require("dseams")`, clone
  [yodaStruct](https://github.com/d-SEAMS/yodaStruct) and run `luarocks make` in the checkout
- R: `remotes::install_github("d-SEAMS/dseams-r")` from [dseams-r](https://github.com/d-SEAMS/dseams-r)
- Julia: `Pkg.add(url="https://github.com/d-SEAMS/Dseams.jl")` from
  [Dseams.jl](https://github.com/d-SEAMS/Dseams.jl)
- Rust: `cargo add dseams@2.14.0` from [dseams-rs](https://github.com/d-SEAMS/dseams-rs)
- Fortran: fpm Git dependency `dseams` from [dseams-fortran](https://github.com/d-SEAMS/dseams-fortran)

[The 2.14.0 notes](reference/release-2.14.md) are the changes
since 2.13.0.
[Migrate from 1.x](user/howto/migrate-1x.md) maps `yodaStruct -c`
onto `seams`.

## Suite stack {#suite-stack}

`libyodaLib` reads a LAMMPS dump or an XYZ file. chemfiles
(`.pdb`, `.gro`, `.dcd`) and the `.con` reader are compiled in when
those libraries are found. The cutoff neighbour list is linkcell when
that library is built, and vesin otherwise. `seams --features` prints
which backends this binary has. The `seams` command, `libseams`, and
`libseams_lang` each link `libyodaLib`. `include/seams.hpp` links
`libseams`. Python and Julia link both libraries. R links
`libseams_lang`. Lua, Rust, and Fortran link `libseams`.

## Documentation structure {#documentation-structure}

The user guide groups the tutorials and the how-tos by
task. The explanation pages are the background for those tasks.

Install the engine and classify one shipped frame.

Tutorials and how-tos, grouped by the analysis you want to run.

The C API, the C++ API, the command line, and the language bindings.

Release notes and how to cite d-SEAMS.

The C++ header reference is
[*api-cpp*](https://docs.dseams.info/api-cpp/).

## Performance {#performance}

The numbers for OVITO, freud, ClavisCage TRACE, SQQ, and PLUMED are on
[Head-to-head](user/comparisons.md).
The sitting medians are on
[the sitting tables](explanation/performance.md).

## Related projects {#related-projects}

**[PydSEAMSlib docs](https://d-seams.github.io/PydSEAMSlib/)**

Python, `import dseams` and `import pydseams`

**[yodaStruct docs](https://d-seams.github.io/yodaStruct/)**

Lua and Fennel, `require("dseams")`

**[linkcell](https://github.com/d-SEAMS/linkcell)**

periodic linked-cell neighbours

**[dseams.info](https://dseams.info)**

project site

## License {#license}

MIT. Cite the 2020 JCIM paper (DOI
[10.1021/acs.jcim.0c00031](https://doi.org/10.1021/acs.jcim.0c00031)).
See [Citation](explanation/citation.md).
