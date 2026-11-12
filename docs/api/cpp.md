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

CI runs the C API test, which is the same binary a C++ caller links:

```bash
meson test -C build seams_c_api --print-errorlogs
```
