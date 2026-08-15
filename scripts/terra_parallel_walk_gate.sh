#!/usr/bin/env bash
# Reconfigure, compile the frame-parallel walk, Catch2, I/O and cages benches.
set -euo pipefail

ROOT="${SEAMS_ROOT:-$HOME/seams-core-merge}"
BUILD="${SEAMS_BUILD:-$ROOT/bcpp}"
PIXI_ENV="${ROOT}/.pixi/envs/default"
DUMP3="${ROOT}/repro/public-ice/nucleation/Ice_nucleation_trajectory/Nucleation_trajectory-3"
LOG="${SEAMS_LOG:-/tmp/seams-parallel-walk-gate.log}"

export PATH="${PIXI_ENV}/bin:${HOME}/.pixi/bin:/usr/bin:/bin"
export LD_LIBRARY_PATH="${BUILD}/src${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"

{
  echo "=== host $(hostname) $(date -Iseconds) ==="
  meson setup --reconfigure "${BUILD}" "${ROOT}"
  meson compile -C "${BUILD}" \
    yodaLib seams test_seams_input bench_lammps_io bench_lammps_walk
  meson test -C "${BUILD}" seams_input seams_read --print-errorlogs

  echo "=== io dump-3 serial ==="
  "${BUILD}/tests/bench_lammps_io" "${DUMP3}" 0 1 1
  echo "=== io dump-3 8 threads ==="
  "${BUILD}/tests/bench_lammps_io" "${DUMP3}" 0 1 8

  echo "=== cages dump-3 first 30 serial ==="
  "${BUILD}/tests/bench_lammps_walk" "${DUMP3}" 30 1 1
  echo "=== cages dump-3 first 30, 8 threads ==="
  "${BUILD}/tests/bench_lammps_walk" "${DUMP3}" 30 1 8

  echo "=== done $(date -Iseconds) ==="
} 2>&1 | tee "${LOG}"
