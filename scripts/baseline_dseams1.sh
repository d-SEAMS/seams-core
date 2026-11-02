#!/usr/bin/env bash
# Build the published d-SEAMS 1 release and time the same cage path as
# tests/bench_cages. The release is tag v1.0.1, commit
# 6e51956f3b9059729f36bcc208f6b6d18c3fec22 (2020-05-21), the 1.x tag
# that carries the JCIM citation 10.1021/acs.jcim.0c00031. Commit
# 1c8efa38 is a 2026 legacy-algorithm surrogate and is not used here.
#
# Paths are relative to this repository unless overridden:
#   SEAMS_ROOT     checkout that contains this script (default: parent)
#   SEAMS_BASELINE worktree of v1.0.1 (default: $SEAMS_ROOT/build/baseline-v1.0.1)
#   SEAMS_TRAJ     LAMMPS dump (default: $SEAMS_ROOT/input/traj/mW_cubic.lammpstrj)
#   SEAMS_REPS     repetitions (default: 1)
#   SEAMS_OUT      JSON path (default: stdout)
set -euo pipefail

SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
ROOT=${SEAMS_ROOT:-$(cd -- "$SCRIPT_DIR/.." && pwd)}
COMMIT=${SEAMS_BASELINE_COMMIT:-6e51956f3b9059729f36bcc208f6b6d18c3fec22}
TREE=${SEAMS_BASELINE:-$ROOT/build/baseline-v1.0.1}
TRAJ=${SEAMS_TRAJ:-$ROOT/input/traj/mW_cubic.lammpstrj}
REPS=${SEAMS_REPS:-1}
FRAME=${SEAMS_FRAME:-1}
BUILD=$TREE/driver-build

if [[ ! -d $ROOT/.git && ! -f $ROOT/.git ]]; then
  echo "SEAMS_ROOT=$ROOT is not a git checkout" >&2
  exit 1
fi
if ! git -C "$ROOT" cat-file -e "$COMMIT^{commit}" 2>/dev/null; then
  echo "commit $COMMIT is not in $ROOT" >&2
  exit 1
fi
GOT=$(git -C "$ROOT" rev-parse "$COMMIT^{commit}")
if [[ $GOT != "$COMMIT" ]]; then
  echo "refusing $GOT; the recorded d-SEAMS 1 commit is $COMMIT" >&2
  exit 1
fi

mkdir -p "$(dirname "$TREE")"
if [[ ! -f $TREE/src/topo_bulk.cpp ]]; then
  git -C "$ROOT" worktree add --detach "$TREE" "$COMMIT"
fi
HEAD=$(git -C "$TREE" rev-parse HEAD)
if [[ $HEAD != "$COMMIT" ]]; then
  echo "worktree $TREE is $HEAD, not $COMMIT" >&2
  exit 1
fi

cxx_ok() {
  printf '%s\n' '#include <iostream>' '#include <boost/multi_array.hpp>' 'int main(){return 0;}' \
    | "$1" -x c++ -fsyntax-only - >/dev/null 2>&1
}
if [[ -n ${CXX:-} ]] && cxx_ok "$CXX"; then
  :
elif cxx_ok g++; then
  CXX=g++
elif [[ -x /usr/bin/g++-13 ]] && cxx_ok /usr/bin/g++-13; then
  CXX=/usr/bin/g++-13
else
  echo "no C++ compiler can see Boost; install libboost-dev" >&2
  exit 1
fi
mkdir -p "$BUILD"
SRCS=(
  mol_sys.cpp
  generic.cpp
  seams_input.cpp
  seams_output.cpp
  neighbours.cpp
  franzblau.cpp
  ring.cpp
  topo_bulk.cpp
  bond.cpp
)
OBJS=()
for src in "${SRCS[@]}"; do
  obj=$BUILD/${src%.cpp}.o
  "$CXX" -std=c++17 -O2 -I"$TREE/src/include/internal" -I"$TREE/src/include/external" \
    -c "$TREE/src/$src" -o "$obj"
  OBJS+=("$obj")
done
echo "baseline compiler: $CXX" >&2
"$CXX" -std=c++17 -O2 -I"$TREE/src/include/internal" -I"$TREE/src/include/external" \
  -o "$BUILD/baseline_dseams1" "$SCRIPT_DIR/baseline_dseams1.cpp" "${OBJS[@]}" \
  -L/usr/lib/x86_64-linux-gnu -lgsl -lgslcblas -lm

SHA=$(sha256sum "$TRAJ" | awk '{print $1}')
JSON=$("$BUILD/baseline_dseams1" "$TRAJ" "$FRAME" "$REPS" 1 "$COMMIT")
# Attach the input checksum without a second parser.
JSON=${JSON/\"input\":/\"input_sha256\": \"$SHA\",$'\n'  \"input\":}
if [[ -n ${SEAMS_OUT:-} ]]; then
  printf '%s\n' "$JSON" > "$SEAMS_OUT"
else
  printf '%s\n' "$JSON"
fi
