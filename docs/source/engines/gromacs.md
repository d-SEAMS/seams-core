# GROMACS

`engines/gromacs/seams_ice.cpp` is a `TrajectoryAnalysisModule`. The binary
is `gmx-seams`. It reads a trajectory with the GROMACS reader. `analyzeFrame`
is the trajectory-analysis callback: it converts nm to Angstrom, calls
`seams_chill_plus`, and records the time. The runner invokes that callback
once per frame.

```bash
gmx-seams -f engines/gromacs/ice.gro -cutoff 3.5
```

`meson test gromacs_mdrun` runs `engines/gromacs/run_mdrun.sh`, and
`meson test docs_engine_gromacs` runs this `-f` / `-cutoff` invocation on
`input/traj/ice1h.lammpstrj`. `run_mdrun.sh`
classifies the trajectory `gmx mdrun` just wrote. Positions are multiplied
by 10 on the way into `seams_chill_plus`. The box conversion is
`gromacs_box_to_lammps` with the same scale.

```{literalinclude} ../../../engines/gromacs/run_mdrun.sh
:language: bash
:lines: 16-21
```

```{literalinclude} ../../../engines/gromacs/seams_ice.cpp
:language: cpp
:start-at: void analyzeFrame
:end-before: void finishAnalysis
```

`ice.gro` is two frames of the 64-atom cubic cell. The tool prints
`overhead_s`, the mean seconds per frame, and for the last frame
`gromacs hexagonal N six-rings M q6 Q`. A topology is not required.
The same step is what an mdrun analysis hook would call. Biasing from inside
mdrun goes through the PLUMED collective variable, which links the same C API.

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

`meson test gromacs_e2e` builds a TIP4P/Ice ice Ic slab with a liquid
layer, runs 200 steps of `gmx mdrun`, classifies the oxygens of every
`.trr` frame with `gmx-seams`, and compares the labels with
`seams chill-plus` on the dump and the dumped oxygens with `gmx dump`.
GROMACS will not read a trajectory from a pipe, so this hook runs after
mdrun. `meson test gromacs_plumed_e2e` instead loads the PLUMED
collective variable into `gmx mdrun -plumed`, classifies inside the MD
loop, and also requires the same labels from `gmx-seams` on the same
steps. With GROMACS 2024 and PLUMED 2.9 that test sets
`GMX_DD_SINGLE_RANK=0`: under single-rank domain decomposition about 5%
of the oxygens reached PLUMED with a neighbouring atom's coordinates.
The test reads the action's timing from `md.log` under the PLUMED patch and
from `PLUMED.OUT` under GROMACS 2025 and later, which link PLUMED without it.

