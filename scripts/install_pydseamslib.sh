#!/bin/sh
# Copyright (c) d-SEAMS developers
# Install pydseamslib from d-SEAMS/PydSEAMSlib, built against this checkout.
#
# The binding's seams-core subproject is this tree, so the docs tests run
# the Python package against the engine under test, before either is
# tagged. PYDSEAMSLIB_REF picks the PydSEAMSlib branch or tag (main).
# Pass --if-missing first to return immediately when pydseamslib is already
# importable in the active Python environment.
# Extra arguments go to pip, e.g. a different interpreter's -m pip.
set -eu
root=$(CDPATH= cd -- "$(dirname "$0")/.." && pwd)
tree=$(git -C "$root" rev-parse HEAD^{tree})
if [ "${1:-}" = "--if-missing" ]; then
  shift
  if SEAMS_ROOT="$root" SEAMS_TREE="$tree" python - <<'PY'
import json
import os
import sys
from pathlib import Path

root = Path(os.environ["SEAMS_ROOT"]).resolve()
expected_tree = os.environ["SEAMS_TREE"]
import dseams  # noqa: F401
import pydseams
if not hasattr(pydseams.Frame, "classify_chill_plus"):
    raise SystemExit("pydseams probe: Frame.classify_chill_plus is missing")
stamp = Path(pydseams.__file__).resolve().parent / ".seams-core-build.json"
if not stamp.is_file():
    raise SystemExit(f"pydseams probe: missing install stamp {stamp}")
info = json.loads(stamp.read_text(encoding="utf-8"))
stamp_root = Path(info.get("root", "")).resolve()
stamp_tree = str(info.get("tree", ""))
if stamp_root != root:
    raise SystemExit(f"pydseams probe: expected root {root}, found {stamp_root}")
if stamp_tree != expected_tree:
    raise SystemExit(f"pydseams probe: expected tree {expected_tree}, found {stamp_tree}")
PY
  then
    exit 0
  fi
fi
ref=${PYDSEAMSLIB_REF:-main}
work=$(mktemp -d)
trap 'rm -rf "$work"' EXIT
git clone --quiet --depth 1 --branch "$ref" https://github.com/d-SEAMS/PydSEAMSlib.git "$work/PydSEAMSlib"
ln -s "$root" "$work/PydSEAMSlib/subprojects/seams-core"
# Meson takes a shared dependency from the top project's wrap, so the
# binding's pins would override the versions this tree asks for.
for wrap in "$root"/subprojects/*.wrap; do
  name=${wrap##*/}
  if [ -f "$work/PydSEAMSlib/subprojects/$name" ]; then
    cp "$wrap" "$work/PydSEAMSlib/subprojects/$name"
  fi
done
cp -R "$root/subprojects/packagefiles/." "$work/PydSEAMSlib/subprojects/packagefiles/"
python -m pip install "$@" "$work/PydSEAMSlib"
SEAMS_ROOT="$root" SEAMS_TREE="$tree" python - <<'PY'
import json
import os
from pathlib import Path

import pydseams

stamp = Path(pydseams.__file__).resolve().parent / ".seams-core-build.json"
stamp.write_text(
    json.dumps({"root": os.environ["SEAMS_ROOT"], "tree": os.environ["SEAMS_TREE"]}, sort_keys=True)
    + "\n",
    encoding="utf-8",
)
PY
