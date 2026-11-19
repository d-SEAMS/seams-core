<!-- https://docs.dseams.info/fronts/gromacs.md -->

# GROMACS tutorial {#gromacs-tutorial}

`gmx-seams` reads a trajectory after `gmx mdrun` and classifies it with `seams_chill_plus`. Coordinates in the trajectory are nanometers. The tool multiplies them by 10 before the call. The lower-triangular GROMACS box stores its rows as the box vectors. The tilts are `box[1][0]`, `box[2][0]` and `box[2][1]` (xy, xz, yz). GROMACS keeps `box[0][1]`, `box[0][2]` and `box[1][2]` at zero.

Build the tool when libgromacs and its trajectory-analysis header are found. The binary is `build/engines/gmx-seams`. `libseams` is `build/src/libseams.so`.

```bash
meson setup build -Dwith_gromacs=auto
meson compile -C build gmx-seams
engines/gromacs/run_mdrun.sh build/engines/gmx-seams build/src
```

`engines/gromacs/run_mdrun.sh` writes a `.tpr` from `topol.top` and `md.mdp`, runs four frozen steps with `gmx mdrun -nsteps 4`, and classifies `ice.xtc`. That trajectory is five frames of the 64-atom cubic cell. The census line is:

```default
gromacs frames 5 ice 320 atoms 320 overhead_s S
```

`ice 320` is the cubic total over those five frames. `overhead_s` is the mean seconds per frame on this machine. `analyzeFrame` is the trajectory-analysis callback that does the conversion and the CHILL+ call.
