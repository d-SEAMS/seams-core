#!/usr/bin/env bash
# Compile the live dump session on the existing meson tree, run Catch2,
# time Niu dump I/O, and walk a short incremental path.
set -euo pipefail

ROOT="${SEAMS_ROOT:-$HOME/seams-core-merge}"
BUILD="${SEAMS_BUILD:-$ROOT/bcpp}"
PIXI_ENV="${ROOT}/.pixi/envs/default"
DUMP_DIR="${ROOT}/repro/public-ice/nucleation/Ice_nucleation_trajectory"
DUMP3="${DUMP_DIR}/Nucleation_trajectory-3"
DUMP2="${DUMP_DIR}/Nucleation_trajectory-2"
LOG="${SEAMS_LOG:-/tmp/seams-dump-session-gate.log}"

export PATH="${PIXI_ENV}/bin:${HOME}/.pixi/bin:/usr/bin:/bin"
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
    test_seams_input \
    test_cage_affiliation \
    bench_lammps_io \
    bench_trajectory_incremental

  echo "=== catch2 ==="
  meson test -C "${BUILD}" seams_input cage_affiliation --print-errorlogs

  echo "=== io dump-3 (O-only x y z, 1551 frames) ==="
  "${BUILD}/tests/bench_lammps_io" "${DUMP3}" 0 1

  echo "=== io dump-2 first 50 (4-site xu yu zu) ==="
  "${BUILD}/tests/bench_lammps_io" "${DUMP2}" 50 1

  echo "=== incremental dump-3 first 30 ==="
  "${BUILD}/tests/bench_trajectory_incremental" "${DUMP3}" 30 1 6 | tee /tmp/niu-inc-30.txt
  echo "=== incremental dump-3 last 3 lines ==="
  tail -5 /tmp/niu-inc-30.txt

  echo "=== done $(date -Iseconds) ==="
} 2>&1 | tee "${LOG}"
