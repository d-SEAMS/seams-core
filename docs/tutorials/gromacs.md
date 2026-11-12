# GROMACS tutorial

d-SEAMS 2.14.0. `gmx-seams` is the trajectory-analysis tool. Build it when
libgromacs is found:

```bash
meson setup build -Dwith_gromacs=auto
meson compile -C build gmx-seams
```

`engines/gromacs/run_mdrun.sh` builds a `.tpr` from `topol.top` and `md.mdp`,
runs four frozen steps with `gmx mdrun`, and classifies `ice.xtc`. CI runs
that script as the `gromacs_mdrun` meson test when the tool was built:

```bash
engines/gromacs/run_mdrun.sh build/gmx-seams build/src
meson test -C build gromacs_mdrun --print-errorlogs
```

The tool prints one line of the form:

```
gromacs frames N ice I atoms A overhead_s S
```

`analyzeFrame` is the trajectory-analysis callback. Coordinates in the
trajectory are nanometres; the tool converts them to Angstrom before
`seams_chill_plus`. The lower-triangular GROMACS box is stored with the
tilts in `box[0][1]`, `box[0][2]` and `box[1][2]` (xy, xz, yz).
