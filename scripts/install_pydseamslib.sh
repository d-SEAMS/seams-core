#!/bin/sh
# Copyright (c) d-SEAMS developers
# Install pydseamslib from d-SEAMS/PydSEAMSlib, built against this checkout.
#
# The binding's seams-core subproject is this tree, so the docs tests run
# the Python package against the engine under test, before either is
# tagged. PYDSEAMSLIB_REF picks the PydSEAMSlib branch or tag (main).
# Extra arguments go to pip, e.g. a different interpreter's -m pip.
set -eu
root=$(CDPATH= cd -- "$(dirname "$0")/.." && pwd)
ref=${PYDSEAMSLIB_REF:-main}
work=$(mktemp -d)
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
