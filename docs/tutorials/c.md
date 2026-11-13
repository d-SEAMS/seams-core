# C tutorial

`seams.h` is the stable API. This program labels one cubic-ice cell and prints the version the library was built with.

```{literalinclude} ../../../docs/examples/ice.c
:language: c
```

```bash
cc -Iinclude -Ibuild/src ice.c -Lbuild/src -lseams -Wl,-rpath,build/src -lm
LD_LIBRARY_PATH=build/src ./a.out
```

`meson test -C build docs_ice_example` runs that program. The binary is `build/src/docs_ice_example`. `tests/test_seams_api.c` is the same call on a 64-atom cell, plus neighbours, q4/q6, and rings. `meson test -C build seams_c_api` runs it. The binary is `build/src/test_seams_api`.


## Structure still

```{figure} ../tutorials/molecules/c/overlay-web.png
:alt: Molecule still for this front
:width: 480px

Ice Ih coloured by published CHILL+, drawn with chemparseplot. `python docs/scripts/render_molecules.py --slug c` writes the 2400 by 1800 PNG.
```

