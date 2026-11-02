#!/usr/bin/env bash
# Reconfigure, compile the frame-parallel walk, Catch2, I/O and cages benches.
set -euo pipefail

SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
ROOT="${SEAMS_ROOT:-$(cd -- "$SCRIPT_DIR/.." && pwd)}"
BUILD="${SEAMS_BUILD:-$ROOT/build}"
PIXI_ENV="${SEAMS_PIXI_ENV:-$ROOT/.pixi/envs/default}"
DUMP_DIR="${SEAMS_NIU_DIR:-$ROOT/data/public-ice/nucleation/Ice_nucleation_trajectory}"
DUMP3="${DUMP_DIR}/Nucleation_trajectory-3"
LOG="${SEAMS_LOG:-$ROOT/logs/seams-parallel-walk-gate.log}"
mkdir -p "$(dirname "$LOG")"

if [[ ! -f "$DUMP3" ]]; then
  python3 "$SCRIPT_DIR/fetch_niu_trajectory.py" --dest "$(dirname "$(dirname "$DUMP_DIR")")"
fi

if [[ -d "${PIXI_ENV}/bin" ]]; then
  export PATH="${PIXI_ENV}/bin:${PATH}"
fi
export LD_LIBRARY_PATH="${BUILD}/src${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"

{
  echo "=== host $(hostname) $(date -Iseconds) ==="
  meson setup --reconfigure "${BUILD}" "${ROOT}"
  meson compile -C "${BUILD}" \
    yodaLib test_seams_input bench_lammps_io bench_lammps_walk
  meson test -C "${BUILD}" seams_input --print-errorlogs

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
