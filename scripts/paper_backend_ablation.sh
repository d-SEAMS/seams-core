#!/usr/bin/env bash
# Same frames as the paper benches. OMP_NUM_THREADS=1 except the OpenMP sweep.
set -euo pipefail
ROOT="${1:-.}"
cd "$ROOT"
export OMP_NUM_THREADS=1
export OMP_PROC_BIND=close
export OMP_PLACES=cores

run_cages() {
  local bin=$1
  local label=$2
  if [[ -x $bin ]]; then
    echo "======== $label ========"
    (cd input && "../$bin" traj/mW_cubic.lammpstrj 1 3)
  else
    echo "======== $label MISSING $bin ========"
  fi
}

run_cages bcpp/tests/bench_cages "default (vesin, closed-form Ylm, Horn)"
run_cages bcpp-sph/tests/bench_cages "sphericart"
run_cages bcpp-ira/tests/bench_cages "ira"
run_cages bcpp-nauty/tests/bench_cages "nauty"

if [[ -x bcpp/tests/bench_overhead ]]; then
  echo "======== vesin overhead N=4096,16000 ========"
  ./bcpp/tests/bench_overhead 4096 16000
fi

if [[ -x bcpp/tests/bench_strong ]]; then
  echo "======== OpenMP steinhardt N=65536 ========"
  for t in 1 2 4; do
    echo "-- OMP_NUM_THREADS=$t --"
    OMP_NUM_THREADS=$t ./bcpp/tests/bench_strong 65536 3
  done
fi
echo DONE
