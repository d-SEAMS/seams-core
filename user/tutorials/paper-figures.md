<!-- https://docs.dseams.info/user/tutorials/paper-figures.md -->

# Paper figures {#paper-figures}

The 2020 JCIM figures are drawn from figshare project 73545.
This page runs the same calls on the supplied example trajectories.
Python is `pip install pydseamslib`; the import is `dseams.userframe`.
Run the block from the repository root.

| figure | trajectory | call | result |
| --- | --- | --- | --- |
| cubic ice, CHILL+ | mW_cubic.lammpstrj frame 1, type 1 | chill_plus() | 4096 cubic |
| cubic ice, cages | the same frame | cages() | n_ic 4096 |
| sI 5^12 cages | =genice_sI.lammpstrj~, type 1 | cages_by_signature("512") | 2 cages, n_atoms 40 |
| mixed frame | =exampleTraj.lammpstrj~, type 2 | cutoff CHILL+ versus cages | interClathrate 9, water 241, cage n_ih 0, n_ic 0, n_water 250 |

`Frame.chill_plus` defaults to the cutoff neighbour list.
On the cubic mW frame that list and `chill_plus(graph="mutual")` both
return 4096 cubic. On the mixed frame the cutoff list is the row above.
The cage score has no hexagonal cage and no double-diamond cage, so
every oxygen is `n_water`.

The nucleation, nanotube, monolayer, and in-plane RDF figures use
trajectories from figshare project 73545.

```python
from dseams import userframe as ds

ice = ds.read("input/traj/mW_cubic.lammpstrj", frame=1, atom_type=1)
print(ice.chill_plus().cubic, ice.chill_plus(graph="mutual").cubic, ice.cages().n_ic)
assert ice.chill_plus().cubic == 4096
assert ice.chill_plus(graph="mutual").cubic == 4096
assert ice.cages().n_ic == 4096

guest = ds.read("input/traj/genice_sI.lammpstrj", atom_type=1)
found = guest.cages_by_signature("512")
print(len(found), found.n_atoms)
assert len(found) == 2 and found.n_atoms == 40

mixed = ds.read("input/traj/exampleTraj.lammpstrj", atom_type=2)
labels = mixed.chill_plus()
score = mixed.cages()
print(labels.interClathrate, labels.water, score.n_ih, score.n_ic, score.n_water)
assert labels.interClathrate == 9 and labels.water == 241
assert score.n_ih == 0 and score.n_ic == 0 and score.n_water == 250
```
