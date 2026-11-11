# CHILL+ on a cubic lattice

The Lua module lives in yodaStruct.
This repository is the C++ engine, and it does not include that module.

Figshare stores the CHILL LAMMPS trajectory of 4096 type-1 mW sites.
The lattice is ice Ic.
CHILL+ labels every site cubic on that lattice.
`nucleation.lammpstrj` is a different figshare deposit.

The 1.x driver is gone.
We use the `seams` command.

## CLI

```bash
seams chill-plus mW_cubic.lammpstrj --cutoff 3.5 --type 1
seams ovito-chill+ mW_cubic.lammpstrj --cutoff 3.5 --type 1
seams chill-plus mW_cubic.lammpstrj --rule examples/ovito-chill-plus.yaml
seams chill-plus mW_cubic.lammpstrj --rule examples/chill-plus.lua
```

`--rule FILE` reads the file as text and takes the first `rule` or
`chill-plus` key.
Comments are skipped.
The file is not executed as Lua or YAML.
A bare `chill+` or `ovito-chill+` is also a `--rule` value.
`ovito-chill+` keeps a shell of four neighbours.
It averages q3 over that whole shell.
It uses open correlation edges and skips the neighbour test.
The default shell count is a bug fix.
The cutoff row is counted before the bond list is cut to four.
Default q3 uses the four nearest neighbours of that four-neighbour shell.

Counts print to stdout.

## Python

`bindings/python/chill_plus.py` is a ctypes wrapper around `libseams_lang`.
It is installed under the data directory `seams/python`.
Set `SEAMS_LANG_LIB` to the built library.
`rule` is a name or the document text.
It does not open a filesystem path, and it is not `pydseams`.

```python
from chill_plus import chill_plus

labels, counts = chill_plus(positions, (lx, ly, lz), rule="chill+")
labels, counts = chill_plus(positions, (lx, ly, lz), rule="ovito-chill+")
```

See the seams CLI manual.
pydseams and yodaStruct are separate pages.

## References

1. Nguyen, A. H., and Molinero, V. (2015).
   Identification of Clathrate Hydrates, Hexagonal Ice, Cubic Ice, and
   Liquid Water in Simulations: the CHILL+ Algorithm.
   The Journal of Physical Chemistry B, 119(29), 9369-9376.
   doi:10.1021/jp510289t
