<!-- https://docs.dseams.info/user/howto/faq.md -->

# FAQ {#faq}

## How do I install the engine? {#how-do-i-install-the-engine}

[Install the engine](../../getting-started/install.md). From a checkout, `pixi run build`
writes `./bbdir/seams`. `nix build` writes `./result/bin/seams`.

## Where is Python? {#where-is-python}

[PydSEAMSlib](https://github.com/d-SEAMS/PydSEAMSlib).
`pip install 'pydseamslib==2.14.0'`.
`import pydseams` is the Frame package.
`import pydseamslib` is that same package.
`import dseams` is the array module.
The compiled extension is `pydseams.yoda`. `pydseams._core` and
`pydseams.cyoda` are the same object.

## Where is Lua / Fennel? {#where-is-lua-fennel}

[yodaStruct](https://github.com/d-SEAMS/yodaStruct).
`require("dseams")`. `require("yoda")` returns that table.
The meson project builds `dseams_core.so`. It does not build a
`yodaStruct` executable. The 1.x command `yodaStruct -c` is
[Migrate from 1.x](migrate-1x.md).

## Why do CHILL+ and cages disagree? {#why-do-chill-and-cages-disagree}

They answer different questions.
[CHILL+ versus cages](../../explanation/chill-plus-vs-cages.md).

## What does `--graph seeded` mean? {#what-does-graph-seeded-mean}

It is the default for `cages` only: mutual k-nearest seeds and the
union of each atom’s k nearest. `k` is 4 unless `-k` is set.
`chill` and `chill-plus` map that name onto the cutoff list.

## Which `--type` do I pass? {#which-type-do-i-pass}

On a single LAMMPS frame, `--type 0` (the default) keeps type 2 when
any atom has it, otherwise type 1.
Mixed TIP4P dumps here, including `input/traj/exampleTraj.lammpstrj`,
use type 2 for oxygen and type 1 for hydrogen.
Single-site dumps use type 1: `ice1h.lammpstrj`,
`interface_f0.lammpstrj`, and `mW_cubic.lammpstrj`.
Any other mixed dump needs the integer of the species you want; see
[Ionic-liquid dumps](ionic-liquid.md).
`density-z` treats `--type 0` as every atom.

## How does a sheared dump box work? {#how-does-a-sheared-dump-box-work}

The dump stores bound spans plus optional `xy`, `xz`, and `yz`.
The edge lengths are recovered from the spans and those tilt factors.
See [Sheared cells](sheared-cells.md).

## Why is chemfiles disabled on `nix build`? {#why-is-chemfiles-disabled-on-nix-build}

chemfiles is not in the flake closure. `seams --features` prints
`chemfiles: disabled`, and `.pdb`, `.gro`, and `.dcd` fail with
`chemfiles is not in this build`. The default pixi lock enables
chemfiles.

## How do I turn colors off? {#how-do-i-turn-colors-off}

Set `NO_COLOR` to a non-empty value. Set `FORCE_COLOR` to a
non-empty value to force color on.

## Does `seams` walk a trajectory? {#does-seams-walk-a-trajectory}

Yes. `--frame N` is the first frame and `--last M` is the last.
Both are 1-based, and `--last` is inclusive. On a LAMMPS dump, a
`--last` past the end of the file stops at the last frame.
`--jobs J` is the OpenMP worker count for that dump. `--jobs 1` is
serial. `.xyz`, `.pdb`, `.gro`, `.dcd`, and `.con` are read one
frame at a time. `seams --features` prints whether OpenMP is on.

## How do I view a classified dump? {#how-do-i-view-a-classified-dump}

`frame.to_solvis()` needs the solvis extra:
`pip install 'pydseamslib[solvis]==2.14.0'`.
The call on ice Ih is on [Python Frame](../../reference/python-frame.md).
The [PydSEAMSlib documentation](https://d-seams.github.io/PydSEAMSlib/)
is the Frame how-to.
`--per-atom FILE` on `chill`, `chill-plus`, `cages`, `fingerprint`,
and `ions` appends a LAMMPS dump with one extra column. OVITO and
VMD open that file.

## How do I cite this? {#how-do-i-cite-this}

[Citation](../../explanation/citation.md). The preferred citation
is the 2020 JCIM paper.
