# C tutorial

Label the eight-atom cubic cell with `seams.h` and print the library version. The program is `docs/examples/ice.c`. It prints `2.14.0 cubic 8`.

```{eval-rst}
.. include:: ../videos/embeds/c.rst
```

```{literalinclude} ../../../docs/examples/ice.c
:language: c
```

`seams.h` includes the generated `seams_version.h`. Meson compiles the program against `libseams` and records the library path in the binary.

```bash
meson setup build
meson compile -C build docs_ice_example
./build/src/docs_ice_example
```

`docs/examples/core_counts.c` reads `input/traj/ice1h.lammpstrj` and `input/traj/genice_sII.lammpstrj`. It returns 0 when Ice Ih is 128 hexagonal waters and 256 six-rings, and sII is 136 clathrate waters, 144 five-rings and 16 six-rings. The lines it prints are:

```
hexagonal 128
rings6 256
clathrate 136
rings5 144
```

```{literalinclude} ../../../docs/examples/core_counts.c
:language: c
```

The function list, status codes and buffer lengths are the {doc}`C API </stable/c>`.

## Ice Ih

```{figure} ../user/tutorials/molecules/c/overlay-web.png
:alt: Ice Ih colored by published CHILL+
:width: 480px

Ice Ih colored by published CHILL+.
```
