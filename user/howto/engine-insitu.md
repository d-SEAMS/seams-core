<!-- https://docs.dseams.info/user/howto/engine-insitu.md -->

# In-situ CHILL+ {#in-situ-chill}

`InSitu` in `adapters/python/dseams_insitu/engine.py` classifies one
frame during a run. `analyse(positions, box)` and
`observe(step, positions, box, ids)` call `seams_frame_chill_rings`.
`observe` also takes `types` and `copy_s`. `analyse_file(path, frame_index)`
reads a LAMMPS dump and calls `seams_frame_chill_plus` and
`seams_frame_rings`. `libseams_lang` is part of the default build.

```bash
meson setup build
meson compile -C build
```

A length-3 box is stored as `0, 0, 0, lx, ly, lz, 0, 0, 0`. Nine LAMMPS
numbers are copied through. Any other length is rejected.

The ice Ih dump is `input/traj/ice1h.lammpstrj`, 128 molecules of type 1.

```bash
seams read input/traj/ice1h.lammpstrj
seams chill-plus input/traj/ice1h.lammpstrj --cutoff 3.5
seams rings input/traj/ice1h.lammpstrj --cutoff 3.5
seams steinhardt input/traj/ice1h.lammpstrj --cutoff 3.5 --l 6
```

`seams read` prints:

```text
nop 128 frame 1 box 15.6457 14.7071 18.073
```

`seams chill-plus` prints `nop 128 hexagonal 128`. `seams rings` prints
`nop 128 rings 256` and then `6 256`, so every primitive ring has size

1. `seams steinhardt --l 6` prints `nop 128 q6 0.625855 q6bar 0.562449`.

Mean local q6 is 0.625855.
