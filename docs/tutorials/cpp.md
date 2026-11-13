# C++ tutorial

d-SEAMS 2.14.0.

```bash
meson setup build -Dwith_tests=true
meson compile -C build test_seams_api docs_ice_example
./build/src/test_seams_api
./build/src/docs_ice_example
```

`test_seams_api` is the `seams_c_api` test. `docs_ice_example` is the eight-atom program in `docs/examples/ice.c`.

```bash
meson test -C build seams_c_api --print-errorlogs
```

A program that classifies one frame links `libseams` and calls
`seams_system_create`, `seams_system_set_frame`, and `seams_chill_plus`.
`seams_set_num_threads(0)` restores the saved `omp_get_max_threads()`.
`type <= 0` unifies mixed types on a copy of the frame.
