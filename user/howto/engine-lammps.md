<!-- https://docs.dseams.info/user/howto/engine-lammps.md -->

# LAMMPS {#lammps}

Build `run_hook` against conda-forge LAMMPS from 2024.8.29 up through
the 2025 releases. Eigen and a C++ compiler need to be on the path.
The exports are the include and library paths of that prefix, including
the AVX directories conda-forge uses.

```bash
pixi init -c conda-forge dseams-lammps
pixi add --manifest-path dseams-lammps/pixi.toml "lammps>=2024.8.29,<2026"
eval "$(pixi shell-hook --manifest-path dseams-lammps/pixi.toml)"
export PKG_CONFIG_PATH="$CONDA_PREFIX/lib/pkgconfig:$CONDA_PREFIX/lib.AVX2_256/pkgconfig:$CONDA_PREFIX/lib.AVX_256/pkgconfig${PKG_CONFIG_PATH:+:$PKG_CONFIG_PATH}"
export LIBRARY_PATH="$CONDA_PREFIX/lib:$CONDA_PREFIX/lib.AVX2_256${LIBRARY_PATH:+:$LIBRARY_PATH}"
export CPATH="$CONDA_PREFIX/include${CPATH:+:$CPATH}"
meson setup build -Dwith_lammps=enabled
meson compile -C build
```

`run_hook` registers a `fix external` callback. LAMMPS passes positions
in angstroms. On a data file written from
`input/traj/ice1h.lammpstrj`, with forces frozen, one callback prints:

```text
lammps hexagonal 128 six-rings 256 q6 0.625855
```

`engines/lammps/in.dseams` is the longer run. The cell is
`engines/lammps/ice.data`, 64 atoms of cubic ice, and the fix line is
`fix ext all external pf/callback 2 2`. With `run 8` the callback runs
5 times. The hook prints `lammps calls 5 stride 2 ice 320` and, for the
last call, `lammps hexagonal 0 six-rings 128 q6 0.628539`. The 320 cubic
labels are 5 calls of 64 atoms. The same coordinates classify as 64
cubic atoms, 128 primitive six-rings, and mean local q6 0.628539.

conda-forge LAMMPS ships `library.h`. This recipe builds `run_hook`.
The compute needs the class headers. The engine page says how to point
Meson at a LAMMPS `src` tree for that plugin.
