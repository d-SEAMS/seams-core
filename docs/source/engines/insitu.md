# In-situ

d-SEAMS 2.14.0. `adapters/python/dseams_insitu/engine.py` is `InSitu`.
`InSitu(library=None, cutoff=3.5, graph="cutoff", max_ring=6, skin=0.0)`
builds one library handle. `graph` is `cutoff` or `mutual`. `skin` above
zero keeps the neighbour list while every oxygen has moved by at most half
that skin. `observe(step, positions, box, ids, types=None, copy_s=0.0)`
classifies one frame through `libseams_lang`. `types` selects the oxygen
indices. There is one OpenMM reporter class, `DseamsReporter`, in
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
meson test -C build insitu_observe docs_engine_insitu --print-errorlogs
```

`observe` stores a length-3 box as `0, 0, 0, lx, ly, lz, 0, 0, 0` and
classifies through `seams_frame_chill_rings`.

```{literalinclude} ../../../engines/insitu_observe.py
:language: python
:lines: 38-40
```

```{literalinclude} ../../../adapters/python/dseams_insitu/engine.py
:language: python
:lines: 199-205
```
