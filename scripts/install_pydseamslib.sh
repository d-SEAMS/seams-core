#!/bin/sh
# Copyright (c) d-SEAMS developers
# Install pydseamslib from d-SEAMS/PydSEAMSlib, built against this checkout.
#
# The binding's seams-core subproject is this tree, so the docs tests run
# the Python package against the engine under test, before either is
# tagged. PYDSEAMSLIB_REF picks the PydSEAMSlib branch or tag (main).
# Pass --if-missing first to return immediately when pydseamslib is already
# installed in the active Python environment for this checkout and HEAD tree.
# Extra arguments go to pip, e.g. a different interpreter's -m pip.
set -eu
root=$(CDPATH= cd -- "$(dirname "$0")/.." && pwd -P)
tree=$(git -C "$root" rev-parse HEAD^{tree})
if [ "${1:-}" = "--if-missing" ]; then
  shift
  if SEAMS_ROOT="$root" SEAMS_TREE="$tree" python -I - <<'PY'
import json
import os
import sysconfig
from pathlib import Path

root = Path(os.environ["SEAMS_ROOT"]).resolve(strict=True)
expected_tree = os.environ["SEAMS_TREE"]
targets = tuple(
    sorted(
        {
            Path(p).resolve()
            for p in (sysconfig.get_path("purelib"), sysconfig.get_path("platlib"))
            if p
        }
    )
)
if not targets:
    raise SystemExit("pydseams probe: could not determine Python site-packages target")


def require_in_target(mod, name):
    mod_file = getattr(mod, "__file__", None)
    if not mod_file:
        raise SystemExit(f"pydseams probe: {name} has no __file__")
    mod_path = Path(mod_file).resolve()
    if not any(mod_path.is_relative_to(target) for target in targets):
        raise SystemExit(
            f"pydseams probe: {name} at {mod_path} is outside install targets {targets}"
        )
    return mod_path


import dseams
import pydseams

require_in_target(dseams, "dseams")
pydseams_path = require_in_target(pydseams, "pydseams")
if not hasattr(pydseams.Frame, "classify_chill_plus"):
    raise SystemExit("pydseams probe: Frame.classify_chill_plus is missing")

stamp = pydseams_path.parent / ".seams-core-build.json"
if stamp.is_symlink() or not stamp.is_file():
    raise SystemExit(f"pydseams probe: missing install stamp {stamp}")
if not any(stamp.resolve().is_relative_to(target) for target in targets):
    raise SystemExit(f"pydseams probe: install stamp {stamp} is outside install targets {targets}")
try:
    info = json.loads(stamp.read_text(encoding="utf-8"))
except (OSError, json.JSONDecodeError) as exc:
    raise SystemExit(f"pydseams probe: unreadable install stamp {stamp}: {exc}")
if not isinstance(info, dict):
    raise SystemExit(f"pydseams probe: malformed install stamp {stamp}")

raw_root = info.get("root")
if not isinstance(raw_root, str) or not raw_root:
    raise SystemExit(f"pydseams probe: missing or invalid stamp root {raw_root!r}")
raw_root_path = Path(raw_root)
if not raw_root_path.is_absolute():
    raise SystemExit(f"pydseams probe: stamp root must be an absolute path, found {raw_root!r}")
stamp_root = raw_root_path.resolve()
if raw_root != str(stamp_root):
    raise SystemExit(
        f"pydseams probe: stamp root must be a canonical realpath, found {raw_root!r} (resolves to {stamp_root})"
    )
if stamp_root != root:
    raise SystemExit(f"pydseams probe: expected root {root}, found {stamp_root}")

stamp_tree = info.get("tree")
if not isinstance(stamp_tree, str) or stamp_tree != expected_tree:
    raise SystemExit(f"pydseams probe: expected tree {expected_tree}, found {stamp_tree!r}")
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
(unset PYTHONPATH; exec python -m pip install "$@" "$work/PydSEAMSlib")
SEAMS_ROOT="$root" SEAMS_TREE="$tree" python -I - <<'PY'
import json
import os
import sysconfig
from pathlib import Path

root = Path(os.environ["SEAMS_ROOT"]).resolve(strict=True)
expected_tree = os.environ["SEAMS_TREE"]
targets = tuple(
    sorted(
        {
            Path(p).resolve()
            for p in (sysconfig.get_path("purelib"), sysconfig.get_path("platlib"))
            if p
        }
    )
)
if not targets:
    raise SystemExit("pydseams stamp: could not determine Python site-packages target")


def require_in_target(mod, name):
    mod_file = getattr(mod, "__file__", None)
    if not mod_file:
        raise SystemExit(f"pydseams stamp: {name} has no __file__")
    mod_path = Path(mod_file).resolve()
    if not any(mod_path.is_relative_to(target) for target in targets):
        raise SystemExit(
            f"pydseams stamp: {name} at {mod_path} is outside install targets {targets}"
        )
    return mod_path


import dseams
import pydseams

require_in_target(dseams, "dseams")
pydseams_path = require_in_target(pydseams, "pydseams")
stamp = pydseams_path.parent / ".seams-core-build.json"
if not any(stamp.parent.resolve().is_relative_to(target) for target in targets):
    raise SystemExit(f"pydseams stamp: {stamp} is outside install targets {targets}")
if stamp.is_symlink():
    stamp.unlink()
stamp.write_text(
    json.dumps({"root": str(root), "tree": expected_tree}, sort_keys=True) + "\n",
    encoding="utf-8",
)
PY
