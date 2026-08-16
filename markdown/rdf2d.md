# In-plane RDF

The Lua scripts live in [yodaStruct](https://github.com/d-SEAMS/yodaStruct).
This repository is the C++ engine.

Figshare:
[In-plane 2D RDF LAMMPS Trajectory](https://figshare.com/articles/In_plane_2D_RDF_LAMMPS_Trajectory/11448711).
fMSI at 320 K. The pair histogram is normalized to an ideal gas.

The `seams` CLI has no RDF command.

## Python

```python
import pydseams as ds

frame = ds.read("dump.lammpstrj")
lines = frame.rdf_2d(output_dir="output/", cutoff=12.0, binwidth=0.05)
```

The engine writes `topoMonolayer/rdf.dat` under `output_dir`.

## Lua

1.x scripts: the RDF helper in
[yodaStruct](https://github.com/d-SEAMS/yodaStruct) `example_lua/`.

## References

1. The in-plane RDF demonstration in the 2020 d-SEAMS paper,
   doi:[10.1021/acs.jcim.0c00031](https://doi.org/10.1021/acs.jcim.0c00031).
