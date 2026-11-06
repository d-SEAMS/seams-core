# One core, one thread, six timed calls after two warmups.
# Writes the ring and neighbour samples and the type 7 tables.
export OMP_NUM_THREADS=1
export OMP_DYNAMIC=FALSE
export OVITO_THREAD_COUNT=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
PY=python3
if [ -x .pixi/envs/default/bin/python ]; then
  PY=.pixi/envs/default/bin/python
fi
STAGES=bbdir/tests/bench_stages
NEIGH=bbdir/tests/bench_neigh_only
if [ ! -x "$STAGES" ]; then
  STAGES=build/tests/bench_stages
fi
if [ ! -x "$NEIGH" ]; then
  NEIGH=build/tests/bench_neigh_only
fi
"$PY" tests/bench_h2h_rings.py \
  --bench-stages-bin "$STAGES" \
  --output-json paper/tex/figures/data/ring-perf/ring-head2head.json \
  --dseams-threads 1 \
  --ovito-threads 1
"$PY" tests/bench_h2h_neighbours.py \
  --bench-bin "$NEIGH" \
  --output-json paper/tex/figures/data/ring-perf/neigh-head2head.json \
  --dseams-threads 1 \
  --freud-threads 1 \
  --ovito-threads 1 \
  --scipy-threads 1
"$PY" paper/tex/figures/ring_perf_type7.py
