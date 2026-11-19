<!-- https://docs.dseams.info/fronts/engines.md -->

# Engines tutorial {#engines-tutorial}

Eight molecules of cubic ice are in `engines/tutorial/ice.data`. Lengths are angstroms.

LAMMPS `run_hook` takes the coordinates LAMMPS already stores in angstroms and calls `seams_chill_plus`. The deck next to the data file is pointed at this cell:

```default
# Eight-molecule cubic ice. The data file is angstroms (units metal).
# fix ext is the in-run hook. The tutorial test classifies ice.data
# with the same seams_chill_plus call the hook makes.
units metal
atom_style atomic
boundary p p p
read_data ice.data
mass 1 16.0
pair_style lj/cut 8.0
pair_coeff * * 0.1 2.75
velocity all set 0.0 0.0 0.0
fix freeze all setforce 0.0 0.0 0.0
fix ext all external pf/callback 1 1
thermo 1
timestep 0.001
run 0
```

The OpenMM plugin reads nanometers. It multiplies positions by 10 before the same call. The tutorial script builds those nanometers from the data file, scales them back, and requires the same eight cubic labels.

```python
def openmm_angstroms(xyz):
    """The plugin reads OpenMM nanometers and multiplies by 10."""
    nanometers = [value / 10.0 for value in xyz]
    return [value * 10.0 for value in nanometers]
```

A passing run prints `engines tutorial ok atoms 8 cubic 8`.
