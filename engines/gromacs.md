<!-- https://docs.dseams.info/engines/gromacs.md -->

# GROMACS {#gromacs}

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

```bash
(
  cd "$work"
  env -u LD_LIBRARY_PATH -u LIBRARY_PATH \
    OMP_NUM_THREADS=1 GMX_MAXBACKUP=-1 \
    gmx grompp -f md.mdp -c ice.gro -p topol.top -o ice.tpr -maxwarn 2
  env -u LD_LIBRARY_PATH -u LIBRARY_PATH \
    OMP_NUM_THREADS=1 GMX_MAXBACKUP=-1 \
    gmx mdrun -s ice.tpr -deffnm ice -nsteps 4 -ntmpi 1
```

```cpp
  void analyzeFrame(int, const t_trxframe &fr, t_pbc *, gmx::TrajectoryAnalysisModuleData *pdata) override {
    if (!fr.bX || !fr.bBox) {
      return;
    }
    const gmx::Selection &sel = pdata->parallelSelection(sel_);
    const int n = sel.posCount();
    if (n <= 0) {
      return;
    }
    std::vector<double> xyz(static_cast<std::size_t>(n) * 3);
    for (int i = 0; i < n; ++i) {
      const auto &x = sel.position(i).x();
      xyz[static_cast<std::size_t>(3 * i)] = static_cast<double>(x[0]) * 10.0;
      xyz[static_cast<std::size_t>(3 * i + 1)] = static_cast<double>(x[1]) * 10.0;
      xyz[static_cast<std::size_t>(3 * i + 2)] = static_cast<double>(x[2]) * 10.0;
    }
    double raw[3][3];
    for (int a = 0; a < 3; ++a) {
      for (int b = 0; b < 3; ++b) {
        raw[a][b] = static_cast<double>(fr.box[a][b]);
      }
    }
    double box9[9];
    gromacs_box_to_lammps(raw, box9, 10.0);
    // frame_census calls seams_chill_plus.
    const auto t0 = std::chrono::steady_clock::now();
    const FrameCensus census = frame_census(xyz.data(), n, box9, cutoffA_);
    const auto t1 = std::chrono::steady_clock::now();
    seconds_ += std::chrono::duration<double>(t1 - t0).count();
    frames_ += 1;
    if (!census.ok) {
      failed_ = true;
      return;
    }
    if (!dump_.empty()) {
      writeDump(fr.bStep ? static_cast<long long>(fr.step) : frames_ - 1, sel, xyz, box9,
                census.labels);
    }
    ice_ += census.cubic;
    atoms_ += n;
    hexagonal_ = census.hexagonal;
    six_ = census.six;
    q6_ = census.q6;
  }
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
