CI builds each engine plugin against the conda-forge engine and runs its `engines-e2e` tests, with LAMMPS on one and two ranks; a skip fails the job.
The PLUMED checks take the `libplumedKernel.so` beside `plumed` when `PLUMED_KERNEL` is unset, and `gromacs_plumed_e2e` reads the action log from `PLUMED.OUT` under GROMACS 2025 and later.
C and Fortran snippets link against the patched chemfiles meson builds when the system one is rejected.
