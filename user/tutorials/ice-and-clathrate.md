<!-- https://docs.dseams.info/user/tutorials/ice-and-clathrate.md -->

# Ice and clathrate labels {#ice-and-clathrate-labels}

You label one frame two ways.
CHILL+ reads a staggered or eclipsed star.
A cage signature counts closed polyhedra of a given ring census.
Install `pydseamslib` with `pip install pydseamslib` and import `from dseams import userframe`, the frame over `libseams_lang`.

## The mixed TIP4P frame {#the-mixed-tip4p-frame}

`input/traj/exampleTraj.lammpstrj` has 750 atoms and 250 oxygens (LAMMPS type 2).

```python
from dseams import userframe as ds

frame = ds.read("input/traj/exampleTraj.lammpstrj", atom_type=2)
print(frame.chill_plus())
print(frame.chill_plus(graph="mutual"))
print(frame.cages())
assert str(frame.chill_plus()) == "IceCounts(interClathrate=9, water=241)"
assert str(frame.chill_plus(graph="mutual")) == "IceCounts(clathrate=1, interClathrate=7, water=242)"
assert str(frame.cages()) == "CageScore(n_ih=0, n_ic=0, n_water=250)"
```

`graph="cutoff"` is the default, the same neighbour list as `seams chill-plus`.
`graph="mutual"` is `seams chill-plus --graph mutual`.
The cage score is the seeded hexagonal and double-diamond membership: no cage is complete, so the water count is 250.
The CLI lines are `nop 250 interClathrate 9 water 241` and `nop 250 clathrate 1 interClathrate 7 water 242`.

## Cubic ice and an sI cage {#cubic-ice-and-an-si-cage}

`input/traj/mW_cubic.lammpstrj` frame 1 is 4096 mW sites.
The file holds 11 frames. `frame=1` selects the first.
Both the cutoff list and the mutual list label every site cubic, and every site sits in a double-diamond cage.

```python
from dseams import userframe as ds

ice = ds.read("input/traj/mW_cubic.lammpstrj", frame=1, atom_type=1)
assert ice.n_atoms == 4096
assert ice.chill_plus().cubic == 4096
assert ice.chill_plus(graph="mutual").cubic == 4096
assert ice.cages().n_ic == 4096
assert ice.cages().n_water == 0
```

`input/traj/genice_sI.lammpstrj` is 46 waters, type 1.
`cages_by_signature("512")` is the $5^{12}$ census: 2 cages and 40 vertices.
`cages_by_signature("51262")` is $5^{12}6^{2}$: 6 cages and 46 vertices.

The default face-edge window is `2.65,3.15` in the coordinate units of the file.
These dumps are angstroms.
`edge="none"` is `--cage-edge none`: every ring the neighbour list already accepted can be a face.
On this sI cell the default and `none` return the same two cages.
The same agreement holds for the sII and sH files in [Classify a clathrate](clathrate.md).
The signature sets the ring depth. An extra `max_ring_size` on this frame is ignored.

```python
from dseams import userframe as ds

frame = ds.read("input/traj/genice_sI.lammpstrj", atom_type=1)
found = frame.cages_by_signature("512")
assert len(found) == 2 and found.n_atoms == 40
found = frame.cages_by_signature("512", edge="none")
assert len(found) == 2 and found.n_atoms == 40
found = frame.cages_by_signature("51262")
assert len(found) == 6 and found.n_atoms == 46
```
