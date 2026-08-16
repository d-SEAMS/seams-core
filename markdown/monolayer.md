# Monolayer square ice

The Lua scripts live in [yodaStruct](https://github.com/d-SEAMS/yodaStruct).
This repository is the C++ engine.

Figshare:
[Monolayer LAMMPS Trajectory](https://figshare.com/articles/Monolayer_LAMMPS_Trajectory/11448741)
(`dump-6-320-310.lammpstrj`, fMSI on cooling 320 K to 310 K).
Four-membered primitive rings.

The `seams` CLI has no monolayer command. Use a front end.

## Python

```python
import pydseams as ds

frame = ds.read("dump-6-320-310.lammpstrj")
frame.monolayer_rings(output_dir="output/", sheet_area=1.0)
```

Pass the sheet area for the system you downloaded.

## Lua

1.x scripts: `example_lua/monolayer/` in
[yodaStruct](https://github.com/d-SEAMS/yodaStruct).

## References

1. The confined-ice topological criterion in the 2020 d-SEAMS paper,
   doi:[10.1021/acs.jcim.0c00031](https://doi.org/10.1021/acs.jcim.0c00031).
