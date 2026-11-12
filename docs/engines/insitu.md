# In-situ

d-SEAMS 2.14.0. `adapters/python/dseams_insitu/engine.py` is `InSitu`.
`observe` classifies one frame through `libseams_lang`. OpenMM's
reporter is the same class: there is one `DseamsReporter`, in
`adapters/python/dseams_insitu/openmm.py` and the driver under
`engines/openmm/reporter.py`.

Build:

```bash
meson setup build -Dwith_tests=true
meson compile -C build
```

CI runs an 8-atom `InSitu.observe` check. It skips when
`libseams_lang.so` or NumPy is missing:

```bash
meson test -C build insitu_observe --print-errorlogs
```
