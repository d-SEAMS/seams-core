# In-situ

d-SEAMS 2.14.0. `adapters/python/dseams_insitu/engine.py` is `InSitu`.
`observe` classifies one frame through `libseams_lang`. There is one
OpenMM reporter class, `DseamsReporter`, in
`adapters/python/dseams_insitu/openmm.py`. It calls `InSitu.observe`.

Build:

```bash
meson setup build -Dwith_tests=true
meson compile -C build
```

CI runs an 8-atom `InSitu.observe` check. The box is three edge lengths,
stored as `0,0,0,lx,ly,lz,0,0,0`. Every label on that cell is cubic (0).
The test skips when `libseams_lang.so` or NumPy is missing, and fails if a
box that is neither 3 nor 9 numbers is accepted:

```bash
meson test -C build insitu_observe --print-errorlogs
```
