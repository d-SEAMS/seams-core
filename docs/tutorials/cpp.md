# C++ tutorial

d-SEAMS 2.14.0.

```bash
c++ -std=c++17 core_counts.cpp -o core-counts $(pkg-config --cflags --libs seams)
./core-counts
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

Ice Ih colored by the published CHILL+ rule.
```
