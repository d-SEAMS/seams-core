<!-- https://docs.dseams.info/engines/lammps.md -->

# LAMMPS {#lammps}

Meson builds `run_hook`, the in-run hook. It links `liblammps` and registers a
`fix external` callback. Every step the callback calls the same C API as
`seams chill-plus`. `engines/lammps/compute_dseams.cpp` is a per-atom style.
When `compute.h`, `atom.h`, and `lammpsplugin.h` are installed, Meson builds
it as `dseamsplugin.so` for `plugin load`. conda-forge LAMMPS installs
`library.h` only; `scripts/fetch_lammps_headers.sh DEST` checks out `src` of
the release that `lmp` reports and prints its path, and `-Dlammps_src_include=`
that path builds the plugin. Otherwise copy the file into a LAMMPS source tree. On several MPI ranks the compute gathers every
atom to rank 0, classifies the whole cell once, and returns each rank its
own labels, so the result does not depend on the domain decomposition.

`compute dseams` takes `cutoff` (angstrom), `graph` (`cutoff` or `mutual`),
`rule` (`published` or `ovito`), `type I` and `mol water`. `type I` scores the
group atoms of LAMMPS type I as water oxygens. `mol water` scores the type
that appears once in each multi-atom molecule and needs an atom style with
molecule IDs. The other group atoms read 6 (unclassified) and atoms outside
the group read -1. A group of several types that sets neither is an error,
so hydrogens and ions stay out of the CHILL+ graph. The deck’s `fix ext` line
is `external pf/callback`.

```default
# mW ice/liquid with compute dseams in situ. The dump carries the same
# coordinates at full precision for the offline CLI.
#   lmp -var plugin dseamsplugin.so -var cells 8 -var steps 200 -var every 10 \
#       -var melt 2000 -var dump insitu.lammpstrj -in in.mw.lmp
include mw_setup.lmp
plugin load ${plugin}
compute ice all dseams cutoff 3.5
dump d all custom ${every} ${dump} id type x y z c_ice
dump_modify d sort id format float %.17g
thermo ${every}
run ${steps}
```

```default
# Cubic ice, metal units (Angstrom). Forces are cleared so the cell does not melt.
# fix ext is the in-run hook. pf/callback 2 2 asks LAMMPS for d-SEAMS every 2 steps.
units metal
atom_style atomic
boundary p p p
read_data ice.data
mass 1 16.0
pair_style lj/cut 8.0
pair_coeff * * 0.1 2.75
velocity all set 0.0 0.0 0.0
fix freeze all setforce 0.0 0.0 0.0
fix ext all external pf/callback 2 2
thermo 2
timestep 0.001
run 8
```

```default
compute ice all dseams
dump 1 all custom 1 dump.chill id type x y z c_ice
run 0
```

```cpp
  double lo[3], hi[3], xy, yz, xz;
  int pflags[3], boxflag;
  lammps_extract_box(hook->lmp, lo, hi, &xy, &yz, &xz, pflags, &boxflag);
  const double box9[9] = {lo[0], lo[1], lo[2], hi[0] - lo[0], hi[1] - lo[1], hi[2] - lo[2],
                          xy, xz, yz};
  // frame_census calls seams_chill_plus.
```
