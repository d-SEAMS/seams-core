The pixi full-suite now fails whenever a test is skipped, so a passing
run confirms that the required analyses and documentation checks all ran.

`-Dwith_lammps=enabled` now requires the LAMMPS class headers needed for
`compute dseams`; Meson errors unless `lammps_src_include` points to
matching sources.
