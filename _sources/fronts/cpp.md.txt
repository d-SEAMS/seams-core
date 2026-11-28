# C++ tutorial

Classify the eight-atom cubic cell with `include/seams.hpp`, then the ice Ih dump. The program is `docs/examples/cpp_api.cpp`. Passed the repository root, the first line is `cpp api cubic 8 q6 0.629` followed by the project version. The second line is:

```
ice Ih atoms 128 hexagonal 128 q6bar 0.562448843702 shell_ih 128 pockets 0
```

```{eval-rst}
.. include:: ../videos/embeds/cpp.rst
```

```{literalinclude} ../../../docs/examples/cpp_api.cpp
:language: cpp
```

`seams.hpp` includes `seams.h`, which includes the generated `seams_version.h`. Meson compiles the program against `libseams`.

```bash
c++ -std=c++20 -Iinclude -Ibuild/src docs/examples/cpp_api.cpp -Lbuild/src -lseams -Wl,-rpath,build/src -o cpp_api
./cpp_api "$PWD"
```

`seams::set_num_threads(0)` restores the OpenMP width captured at the first call, `omp_get_max_threads()`. A type of 0 or less on a frame whose atoms do not all share one positive type classifies a copy in which those atoms are given one type. That copy is `unified_copy` in `src/seams_api.cpp`.

`docs/examples/core_counts.cpp` checks the same Ice Ih and sII census. It prints `cpp hexagonal 128 rings6 256 clathrate 136 rings5 144` and returns 0 only when sII also has 16 six-rings.

```{literalinclude} ../../../docs/examples/core_counts.cpp
:language: cpp
```

The class list is the {doc}`C++ API </stable/cpp>`.

## Ice Ih

```{figure} ../user/tutorials/molecules/cpp/overlay-web.png
:alt: Ice Ih colored by published CHILL+
:width: 480px

Ice Ih colored by published CHILL+.
```
