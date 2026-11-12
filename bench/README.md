# Benchmarks

`seams_bench` times four analyses on one water/ice frame:

| Case | Work |
| --- | --- |
| `neigh` | Cutoff neighbour list (`nneigh::neighListO`) |
| `q3`, `q6` | Steinhardt bond-order parameters, local and neighbour-averaged |
| `rings` | Primitive (Franzblau) rings up to length 6 |
| `ice` | Topological ice identification (`findHC` then `findDDC`) |

The default input is a periodic ice Ic lattice: 8×8×8 conventional
diamond cells, 4096 molecules, lattice constant 6.37 Å, neighbour
cutoff 3.5 Å. That is the same molecule count as
`input/traj/mW_cubic.lammpstrj`. Pass `--traj` to time a dump instead.
Setup (reading the frame, and any analysis a later case consumes) is
outside the timed region. Each printed line is one full run.

Build from a Release tree (the executable lands in `bbdir-rel/bench/`):

```bash
meson setup bbdir-rel --buildtype=release -Ddebug=true -Dwith_tests=true
meson compile -C bbdir-rel seams_bench
OMP_NUM_THREADS=1 ./bbdir-rel/bench/seams_bench --case all --reps 1
```

`bench/measure.py` compares two binaries. It runs at least 20 rounds,
alternating which binary goes first, then a same-binary A/A band and a
bootstrap 95% interval of the median ratio. Outputs go under
`bench/local/` (gitignored), along with `*.bench.csv`, `*.bench.json`,
and `*.bench.log`.

```bash
python3 bench/measure.py \
  --base bbdir-base/bench/seams_bench \
  --patch bbdir-rel/bench/seams_bench \
  --rounds 20
```

`--lib-path` adds directories to `LD_LIBRARY_PATH` when a binary
was copied out of the build tree and still needs the build's
shared libraries. The timed region is one full analysis after one
warmup, on one thread.
