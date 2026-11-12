# GROMACS

`engines/gromacs/seams_ice.cpp` is a `TrajectoryAnalysisModule`. `gmx-seams`
reads a trajectory with the GROMACS reader. `analyzeFrame` converts nm to
Angstrom, calls `seams_chill_plus`, and records the time. That method is the
in-run hook: it is invoked once per frame by the analysis runner.

```bash
gmx-seams -f engines/gromacs/ice.gro -cutoff 3.5
```

`ice.gro` is two frames of the 64-atom cubic cell. The tool prints
`overhead_s`, the mean seconds per frame. A topology is not required.
The same step is what an mdrun analysis hook would call. Biasing from inside
mdrun goes through the PLUMED collective variable, which links the same C API.
