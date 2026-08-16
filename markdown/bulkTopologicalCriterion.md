# Bulk HC / DDC cages

The Lua scripts live in [yodaStruct](https://github.com/d-SEAMS/yodaStruct).
This repository is the C++ engine.

Figshare:
[Nucleation LAMMPS Trajectory](https://figshare.com/articles/Nucleation_LAMMPS_Trajectory/11448702)
(4096 mW waters after crystallization). Hexagonal cages (HC),
double-diamond cages (DDC), and mixed rings sit on the ice cluster.

`seams cages` replaces the 1.x `bulkTopologicalCriterion` Lua script.

## CLI

```bash
seams cages nucleation.lammpstrj --type 2
```

Add `--graph knn` / `knn-union` / `seeded` when the bond graph is not
the cutoff list.

## Python

```python
import pydseams as ds

frame = ds.read("nucleation.lammpstrj")
print(frame.cages())
```

## Lua

```lua
local dseams = require("dseams")
local cloud = dseams.read("nucleation.lammpstrj", {type = 2})
print(dseams.cages(cloud))
```

Books: [seams CLI](https://docs.dseams.info),
[pydseams](https://d-seams.github.io/PydSEAMSlib/).

## References

1. Haji-Akbari, A.; Debenedetti, P. G. *Proc. Natl. Acad. Sci.* **2015**,
   *112*, 10582. doi:[10.1073/pnas.1509267112](https://doi.org/10.1073/pnas.1509267112)
