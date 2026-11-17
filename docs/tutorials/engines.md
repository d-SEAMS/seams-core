# Engines tutorial

LAMMPS `run_hook` takes the coordinates LAMMPS already stores in angstroms and calls `seams_chill_plus`. The deck next to the data file is the same hook the in-tree test uses, pointed at this cell:

```{literalinclude} ../../../engines/tutorial/in.dseams
```

The OpenMM plugin reads nanometers. It multiplies positions by 10 before the same call. The tutorial test builds those nanometers from the data file, scales them back, and requires the same eight cubic labels.

```{literalinclude} ../../../engines/tutorial_walk.py
:language: python
:pyobject: openmm_angstroms
```

A passing run prints `engines tutorial ok atoms 8 cubic 8`.
