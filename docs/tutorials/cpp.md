# C++ tutorial

d-SEAMS 2.14.0.

```bash
meson setup build -Dwith_tests=true
meson compile -C build test_seams_api docs_ice_example
./build/src/test_seams_api
./build/src/docs_ice_example
```

```{literalinclude} ../../../docs/examples/core_counts.cpp
:language: cpp
```

A program that classifies one frame links `libseams` and calls
`seams_system_create`, `seams_system_set_frame`, and `seams_chill_plus`.
`seams_set_num_threads(0)` restores the saved `omp_get_max_threads()`.
`type <= 0` unifies mixed types on a copy of the frame.

## Structure still

```{figure} ../tutorials/molecules/cpp/overlay-web.png
:alt: Molecule still for this front
:width: 480px

Ice Ih colored by published CHILL+, drawn with chemparseplot. `python docs/scripts/render_molecules.py --slug cpp` writes the 2400 by 1800 PNG.
```
