# GROMACS

`engines/gromacs/seams_ice.cpp` is a `TrajectoryAnalysisModule`. The binary
is `gmx-seams`. It reads a trajectory with the GROMACS reader. `analyzeFrame`
is the trajectory-analysis callback: it converts nm to Angstrom, calls
`seams_chill_plus`, and records the time. The runner invokes that callback
once per frame.

```bash
gmx-seams -f engines/gromacs/ice.gro -cutoff 3.5
```

`run_mdrun.sh` classifies the trajectory `gmx mdrun` wrote. Positions are
multiplied by 10 on the way into `seams_chill_plus`. The box conversion is
`gromacs_box_to_lammps` with the same scale.

```{literalinclude} ../../../engines/gromacs/run_mdrun.sh
:language: bash
:lines: 21-28
```

```{literalinclude} ../../../engines/gromacs/seams_ice.cpp
:language: cpp
:start-at: void analyzeFrame
:end-before: void finishAnalysis
```

`ice.gro` is two frames of the 64-atom cubic cell. The tool prints
`overhead_s`, the mean seconds per frame, and for the last frame
`gromacs hexagonal N six-rings M q6 Q`. `overhead_s` is machine-dependent.
A topology is not required.
`gmx-seams` runs after `gmx mdrun` and reads the trajectory that run wrote.
Biasing from inside mdrun goes through the
[d-SEAMS PLUMED plugin v0.4.1](plumed.md), which links the same C API.

`-select` picks the sites to classify; the default is every atom. A
multi-site water model needs its oxygens only, which takes a run input
file: `gmx-seams -s topol.tpr -f traj.trr -select 'name OW'`. `-dump file`
appends a LAMMPS dump frame per analysed frame with the coordinates
handed to d-SEAMS, in Angstrom, and the CHILL+ label. Molecules are not
made whole, so those coordinates are the ones in the trajectory.

GROMACS 2023 and later install neither the trajectory-analysis headers
nor the generated `libgromacs_export.h` they include. Point the
`gromacs_legacy_include` option at `api/legacy/include` of the GROMACS
source that matches the installed `libgromacs`; `engines/gromacs/compat`
supplies the export header. `gmx-seams` builds as C++17, as GROMACS does.
