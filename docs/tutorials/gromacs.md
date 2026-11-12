# GROMACS tutorial

`seams_ice` is a trajectory analysis tool. `engines/gromacs/run_mdrun.sh` builds a `.tpr` from `topol.top` and `md.mdp`, runs four frozen steps with `gmx mdrun`, and classifies `ice.xtc`.

```bash
meson compile -C build seams_ice
engines/gromacs/run_mdrun.sh build/engines/seams_ice build/src
```

The tool prints one line:

```
gromacs frames 5 ice 320 atoms 320 overhead_s ...
```

Five frames, including the start, are entirely cubic. `analyzeFrame` is the per-frame hook. Coordinates in the trajectory are nanometres; the tool converts them to Angstrom before `seams_chill_plus`.
