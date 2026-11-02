#!/usr/bin/env bash
# Compile the live dump session on the existing meson tree, run Catch2,
# time Niu dump I/O, and walk a short incremental path.
set -euo pipefail

SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
ROOT="${SEAMS_ROOT:-$(cd -- "$SCRIPT_DIR/.." && pwd)}"
BUILD="${SEAMS_BUILD:-$ROOT/build}"
PIXI_ENV="${SEAMS_PIXI_ENV:-$ROOT/.pixi/envs/default}"
DUMP_DIR="${SEAMS_NIU_DIR:-$ROOT/data/public-ice/nucleation/Ice_nucleation_trajectory}"
DUMP3="${DUMP_DIR}/Nucleation_trajectory-3"
DUMP2="${DUMP_DIR}/Nucleation_trajectory-2"
LOG="${SEAMS_LOG:-$ROOT/logs/seams-dump-session-gate.log}"
mkdir -p "$(dirname "$LOG")" "$DUMP_DIR"

if [[ ! -f "$DUMP3" || ! -f "$DUMP2" ]]; then
  python3 "$SCRIPT_DIR/fetch_niu_trajectory.py" --dest "$(dirname "$(dirname "$DUMP_DIR")")"
fi

if [[ -d "${PIXI_ENV}/bin" ]]; then
  export PATH="${PIXI_ENV}/bin:${PATH}"
fi
export LD_LIBRARY_PATH="${BUILD}/src${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"

{
  echo "=== host $(hostname) $(date -Iseconds) ==="
  echo "root=${ROOT}"
  echo "build=${BUILD}"
  test -f "${ROOT}/src/seams_input.cpp"
  test -f "${BUILD}/build.ninja"
  test -x "${PIXI_ENV}/bin/meson" || test -x /usr/bin/meson

  echo "=== reconfigure ==="
  meson setup --reconfigure "${BUILD}" "${ROOT}"

  echo "=== compile ==="
  meson compile -C "${BUILD}" \
    yodaLib \
    seams \
    test_seams_input \
    test_cage_affiliation \
    bench_lammps_io \
    bench_lammps_walk \
    bench_trajectory_incremental

  echo "=== catch2 ==="
  meson test -C "${BUILD}" seams_input cage_affiliation --print-errorlogs

  echo "=== io dump-3 (O-only x y z, 1551 frames) ==="
  "${BUILD}/tests/bench_lammps_io" "${DUMP3}" 0 1

  echo "=== io dump-2 first 50 (4-site xu yu zu) ==="
  "${BUILD}/tests/bench_lammps_io" "${DUMP2}" 50 1

  echo "=== incremental dump-3 first 30 ==="
  INC_LOG="${SEAMS_INC_LOG:-$ROOT/logs/niu-inc-30.txt}"
  mkdir -p "$(dirname "$INC_LOG")"
  "${BUILD}/tests/bench_trajectory_incremental" "${DUMP3}" 30 1 6 | tee "$INC_LOG"
  echo "=== incremental dump-3 last 3 lines ==="
  tail -5 "$INC_LOG"

  echo "=== done $(date -Iseconds) ==="
} 2>&1 | tee "${LOG}"
