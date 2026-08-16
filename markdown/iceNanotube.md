# Ice nanotube prisms

The Lua scripts live in [yodaStruct](https://github.com/d-SEAMS/yodaStruct).
This repository is the C++ engine.

Figshare:
[Quasi-1D Nanotube LAMMPS Trajectory](https://figshare.com/articles/Quasi-1D_Nanotube_LAMMPS_Trajectory/11448768)
(`dump-240-square.lammpstrj`, 750 TIP4P/2005 waters in a smooth (13,0)
tube at 240 K). Tetragonal prism blocks.

The `seams` CLI has no prism command. Use a front end.

## Python

```python
import pydseams as ds

frame = ds.read("dump-240-square.lammpstrj")
frame.find_prisms(output_dir="output/")
```

## Lua

The 1.x scripts live under `example_lua/iceNanotube/` in
[yodaStruct](https://github.com/d-SEAMS/yodaStruct). The library call
is `require("dseams")`, not `yodaStruct -c`.

## References

1. The confined-ice topological criterion in the 2020 d-SEAMS paper,
   doi:[10.1021/acs.jcim.0c00031](https://doi.org/10.1021/acs.jcim.0c00031).
