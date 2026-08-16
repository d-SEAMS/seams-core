# CHILL+ on a cubic lattice

The Lua scripts live in [yodaStruct](https://github.com/d-SEAMS/yodaStruct).
This repository is the C++ engine.

Figshare:
[CHILL LAMMPS Trajectory](https://figshare.com/articles/CHILL_LAMMPS_Trajectory/11448720)
(`mW_cubic.lammpstrj`, 4096 type-1 mW sites of ice Ic). CHILL+
labels every site cubic on that lattice. `nucleation.lammpstrj` is
a different figshare deposit.

The 1.x driver (`yodaStruct -c`, `conf.yaml`, `lua_inputs/`) is gone.
Use `seams`, `pydseams`, or `require("dseams")`.

## CLI

```bash
seams chill-plus mW_cubic.lammpstrj --cutoff 3.5 --type 1
```

Counts print to stdout (`cubic`, `hexagonal`, `water`,
`interfacial`, `clathrate`, `interClathrate`, `reCubic`, `reHex`,
`unclassified`).

## Python

```python
import pydseams as ds

frame = ds.read("mW_cubic.lammpstrj")
print(frame.chill_plus())
```

## Lua

```lua
local dseams = require("dseams")
local cloud = dseams.read("mW_cubic.lammpstrj", {type = 1})
print(dseams.chill_plus(cloud, {cutoff = 3.5}))
```

Books: [seams CLI](https://docs.dseams.info),
[pydseams](https://d-seams.github.io/PydSEAMSlib/),
[dseams](https://d-seams.github.io/yodaStruct/).

## References

1. Nguyen, A. H., & Molinero, V. (2014). Identification of Clathrate Hydrates, Hexagonal Ice, Cubic Ice, and Liquid Water in Simulations: the CHILL+ Algorithm. The Journal of Physical Chemistry B, 119(29), 9369-9376. doi:10.1021/jp510289t
