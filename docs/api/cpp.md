# C++

d-SEAMS 2.14.0. C++ callers include `seams.h` and link `libseams`. The
header is C linkage. Coordinates are Angstrom. A cell is nine LAMMPS
numbers: `ox, oy, oz, lx, ly, lz, xy, xz, yz`.

`libseams_lang` is the wider C ABI used by the language fronts
(`seams_frame_*`). It is not required for a program that only needs
`seams.h`.

`docs/examples/minimal_cpp.cpp` reads the
GenIce sI cell with `seams_frame_read`, then calls
`seams_frame_neighbours`, `seams_frame_steinhardt`,
`seams_frame_chill_plus`, and `seams_frame_rings`.
`meson test docs_examples` runs it.

```{literalinclude} ../../../docs/examples/minimal_cpp.cpp
:language: cpp
```

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
`meson test docs_api_cpp` compiles and runs the eight-atom program below.
It calls `seams_set_num_threads`, `seams_chill_plus`, and `seams_q4q6`.

```{literalinclude} ../../../docs/examples/ice.cpp
:language: cpp
```

```{literalinclude} ../../../bindings/cxx/test_lang.cpp
:language: cpp
:lines: 33-40,53-75,108-120,141-159,184-208,301-357,449-453,495-502,1058-1075,1119-1192
```
