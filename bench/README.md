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

The full trajectory is in the public Hugging Face dataset
`d-SEAMS/seams-core-inputs`. Download and verify its SHA-256 digest
with `input/fetch.py`:

```bash
python3 input/fetch.py mW_cubic.lammpstrj
```

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

## C API thread scaling

`capi_thread_scale` calls `seams_cutoff_csr` and `seams_ql_cutoff` on
24×24×24 conventional cells of ice Ic (110592 oxygens, cutoff 3.5 Å, q6)
at 1, 2, 4 and 8 OpenMP threads. Each width must return the one-thread
bytes. `seams_c_last_stages` splits every call into point-cloud
conversion, cutoff-list build, q_l and copy-out; `other_ms` is argument
checks and freeing the temporaries. `capi_thread_scale.txt` here is the
median of 21 calls per width from a release tree:

```bash
./build/capi_thread_scale bench/capi_thread_scale.txt 24 21
```

At 8 threads the cutoff-list build is the larger part of both calls
(about 74% of `seams_cutoff_csr`, 60% of `seams_ql_cutoff`) and speeds up
only 2.2× over one thread. The q_l stage speeds up 5.4× and is under a
fifth of the time. The point-cloud conversion is serial, so it stays
near 0.8 ms at every width and becomes 15–18% of the 8-thread time.
