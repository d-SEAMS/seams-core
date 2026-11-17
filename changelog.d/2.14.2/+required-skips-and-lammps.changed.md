Any skipped test now fails the pixi full-suite CI job. The suite
requires `sqq_sh_triclinic`, the six `docs_paper_*_full` checks, and
`docs_render_molecules` once it installs their dependencies and
enables `paper_full`.

`-Dwith_lammps=enabled` now requires the LAMMPS class headers needed for `compute dseams`; Meson errors unless `lammps_src_include` points to matching sources. CI jobs that enable LAMMPS now provide that include path.
