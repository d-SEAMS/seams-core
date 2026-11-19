<!-- https://docs.dseams.info/engines/plumed.md -->

# PLUMED {#plumed}

The [d-SEAMS PLUMED plugin egg](https://www.plumed-nest.org/eggs/26/018/)
(plumID:26.018) uses v0.4.1 of the
[dseams-plumed repository](https://github.com/d-SEAMS/dseams-plumed).
The build below uses PLUMED 2.9. A conda-forge LAMMPS package links
that release, so PLUMED 2.10 stays in its own environment.

The thread count comes from `PLUMED_NUM_THREADS`. Set that variable
before PLUMED starts. `OMP_NUM_THREADS` does not change the PLUMED team.
The standalone PLUMED action (`dseams-plumed` v0.4.1) is indexed in
[PLUMED-NEST egg 26.018](https://www.plumed-nest.org/eggs/26/018/).

`SEAMS` is a collective variable with five components.

- `seams.q6` is the mean local q6. Its derivatives are a central difference, so
  `RESTRAINT` can bias it. A `PRINT` leaves that derivative array empty, and
  the difference runs only after a `RESTRAINT` has allocated it.
- `seams.ice` is the fraction of atoms published CHILL+ calls cubic. The label
  is discrete, so the component has no derivatives and is not biasable.
- `seams.ring` is the fraction of atoms on a primitive six-ring. It is not
  biasable either: a `RESTRAINT` on `seams.ice` or `seams.ring` is rejected
  because those components carry no derivatives.
- `seams.hexagonal` is the CHILL+ hexagonal count.
- `seams.six` is the primitive six-ring count.

`MODE=q6` computes only the mean local q6 and `MODE=q6bar` only the
mean Lechner-Dellago q6-bar, which has no derivatives. Neither runs CHILL+
or the ring search, and neither has the other four components. `DUMP=file`
appends one LAMMPS dump frame per call with the coordinates handed to
d-SEAMS, in angstroms, the per-atom q6 or q6-bar, and under `MODE=all` the
CHILL+ label. `seams chill-plus` reads that dump back.

```default
LOAD FILE=./libSeamsPlumed.so
seams: SEAMS ATOMS=1-64 CUTOFF=3.5
RESTRAINT ARG=seams.q6 AT=0.5 KAPPA=10.0
PRINT ARG=seams.q6,seams.ice,seams.ring FILE=colvar STRIDE=1
```

```cpp
bool SeamsCV::score(const std::vector<double> &xyz, const double box9[9], Score &out, bool keep) {
  const int n = static_cast<int>(xyz.size() / 3);
  seams_system *sys = seams_system_create();
  if (sys == nullptr || seams_system_set_frame(sys, xyz.data(), n, box9, nullptr, nullptr) != SEAMS_OK) {
    seams_system_destroy(sys);
    return false;
  }
  peratom_.resize(static_cast<std::size_t>(n));
  bool ok = true;
  if (mode_ == Mode::q6bar) {
    double mean = 0.0;
    ok = seams_q6bar(sys, cutoff_, 0, 0.36, peratom_.data(), peratom_.size(), &mean, nullptr,
                     nullptr) == SEAMS_OK;
    double sum = 0.0;
    for (double value : peratom_) sum += value;
    out.q6 = sum / static_cast<double>(n);
  } else {
    ok = seams_q4q6(sys, cutoff_, 0, nullptr, peratom_.data()) == SEAMS_OK;
    double sum = 0.0;
    for (double value : peratom_) sum += value;
    out.q6 = sum / static_cast<double>(n);
  }
  if (ok && mode_ == Mode::all && keep) {
    labels_.resize(static_cast<std::size_t>(n));
    int counts[9] = {};
    seams_rings *rings = nullptr;
    int hist[7] = {};
    ok = seams_chill_plus(sys, cutoff_, 0, SEAMS_GRAPH_CUTOFF, SEAMS_CHILL_PUBLISHED,
                          labels_.data(), counts) == SEAMS_OK &&
         seams_rings_find(sys, cutoff_, 0, SEAMS_GRAPH_CUTOFF, 6, &rings) == SEAMS_OK &&
         seams_rings_fraction(rings, 6, &out.ring) == SEAMS_OK &&
         seams_rings_histogram(rings, hist, 7) == SEAMS_OK;
    seams_rings_destroy(rings);
    out.ice = static_cast<double>(counts[SEAMS_LABEL_CUBIC]) / static_cast<double>(n);
    out.hexagonal = counts[SEAMS_LABEL_HEXAGONAL];
    out.six = hist[6];
  }
  seams_system_destroy(sys);
  return ok;
}
```
