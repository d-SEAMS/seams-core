# GROMACS

`engines/gromacs/seams_ice.cpp` is a `TrajectoryAnalysisModule`. The binary
is `gmx-seams`. It reads a trajectory with the GROMACS reader. `analyzeFrame`
is the trajectory-analysis callback: it converts nm to Angstrom, calls
`seams_chill_plus`, and records the time. The runner invokes that callback
once per frame.

```bash
gmx-seams -f engines/gromacs/ice.gro -cutoff 3.5
```

`meson test gromacs_mdrun` runs `engines/gromacs/run_mdrun.sh`, which
classifies the trajectory `gmx mdrun` just wrote. Positions are multiplied
by 10 on the way into `seams_chill_plus`. The box conversion is
`gromacs_box_to_lammps` with the same scale.

```{literalinclude} ../../../engines/gromacs/run_mdrun.sh
:language: bash
:lines: 16-21
```

```{literalinclude} ../../../engines/gromacs/seams_ice.cpp
:language: cpp
:lines: 44-62
```

`ice.gro` is two frames of the 64-atom cubic cell. The tool prints
`overhead_s`, the mean seconds per frame. A topology is not required.
The same step is what an mdrun analysis hook would call. Biasing from inside
mdrun goes through the PLUMED collective variable, which links the same C API.
