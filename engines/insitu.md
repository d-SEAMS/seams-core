<!-- https://docs.dseams.info/engines/insitu.md -->

# In-situ {#in-situ}

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

`observe` stores a length-3 box as `0, 0, 0, lx, ly, lz, 0, 0, 0` and
classifies through `seams_frame_chill_rings`.

```python
box = np.array([a, a, a], dtype=np.float64)
engine = InSitu()
engine.observe(0, pos, box, ids=np.arange(8))
```

```python
    def _classify(self):
        self._fail(
            self.lib.seams_frame_chill_rings(
                ctypes.c_void_p(self._frame),
                ctypes.c_double(self.cutoff),
                ctypes.c_int(1),
                ctypes.c_char_p(self.graph.encode()),
                ctypes.c_int(self.max_ring),
```
