# C tutorial

`seams.h` is the stable API. This program labels one cubic-ice cell and prints the version the library was built with.

```{literalinclude} ../../../docs/examples/core_counts.c
:language: c
```

```{literalinclude} ../../../docs/examples/ice.c
:language: c
```

```bash
cc -Iinclude -Ibuild/src ice.c -Lbuild/src -lseams -Wl,-rpath,build/src -lm
LD_LIBRARY_PATH=build/src ./a.out
```

## Structure still

```{figure} ../tutorials/molecules/c/overlay-web.png
:alt: Molecule still for this front
:width: 480px

Ice Ih colored by published CHILL+, drawn with chemparseplot. `python docs/scripts/render_molecules.py --slug c` writes the 2400 by 1800 PNG.
```
