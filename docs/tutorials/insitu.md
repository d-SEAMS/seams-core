# In-situ tutorial

d-SEAMS 2.14.0. Eight cubic-ice atoms, one observation:

```bash
meson setup build -Dwith_tests=true
meson compile -C build
python3 engines/insitu_observe.py build
```

That script is the `insitu_observe` meson test CI runs. It exits 77
when NumPy or `libseams_lang.so` is absent, and otherwise checks that
`InSitu.observe` returns eight labels.

```bash
meson test -C build insitu_observe --print-errorlogs
```
