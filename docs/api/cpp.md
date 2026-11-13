# C++

d-SEAMS 2.14.0. C++ callers include `seams.h` and link `libseams`. The
header is C linkage. Coordinates are Angstrom. A cell is nine LAMMPS
numbers: `ox, oy, oz, lx, ly, lz, xy, xz, yz`.

`libseams_lang` is the wider C ABI used by the language fronts
(`seams_frame_*`). It is not required for a program that only needs
`seams.h`.

Build:

```bash
meson setup build -Dwith_tests=true
meson compile -C build
```

CI runs the C API test, which is the same binary a C++ caller links.
The executable is `build/src/test_seams_api`.

```bash
meson test -C build seams_c_api --print-errorlogs
meson test -C build seams_lang_fixtures --print-errorlogs
```

`bindings/cxx/test_lang.cpp` is the `seams_lang_fixtures` test. It calls
every function in `seams_lang.h` on the in-tree trajectories.

```{literalinclude} ../../../bindings/cxx/test_lang.cpp
:language: cpp
:lines: 33-40,53-75,108-120,141-159,184-208,301-357,449-453,495-502,1058-1075,1119-1192
```
