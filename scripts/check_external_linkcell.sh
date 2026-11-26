#!/bin/sh
# Copyright (c) d-SEAMS developers
# Install the pinned linkcell and build exception_typeinfo against that
# prefix. The executable and the hidden-visibility probe must list
# liblinkcell_eh in their own NEEDED entries, then the typeinfo test runs
# with no loader path overriding the prefix.
set -eu

root=$(CDPATH= cd -- "$(dirname "$0")/.." && pwd)
cd "$root"
base=${1:-/tmp/sc-external-linkcell}
mkdir -p "$base"
prefix=$base/prefix
lc_build=$base/linkcell-build
seams_build=$base/seams-build
lc_src=$root/subprojects/linkcell
pin=$(sed -n 's/^revision = //p' "$root/subprojects/linkcell.wrap")
got=""
if git -C "$lc_src" rev-parse HEAD > "$base/linkcell-rev.txt" 2>/dev/null; then
  got=$(cat "$base/linkcell-rev.txt")
fi
if [ ! -f "$lc_src/meson.build" ] || [ "$got" != "$pin" ]; then
  rm -rf "$lc_src"
  meson subprojects download linkcell
  git -C "$lc_src" rev-parse HEAD > "$base/linkcell-rev.txt"
  got=$(cat "$base/linkcell-rev.txt")
fi
if [ "$got" != "$pin" ]; then
  echo "linkcell checkout $got is not the wrap pin $pin" >&2
  exit 1
fi
if [ ! -f "$root/subprojects/rang-3.2/meson.build" ]; then
  meson subprojects download rang
fi
if [ ! -f "$root/subprojects/cxxopts-3.2.0/meson.build" ]; then
  meson subprojects download cxxopts
fi
if [ ! -f "$root/subprojects/Catch2-3.14.0/meson.build" ]; then
  meson subprojects download catch2
fi

meson setup "$lc_build" "$lc_src" \
  --prefix="$prefix" --libdir=lib --buildtype=debug \
  -Dwith_tests=false -Dwith_examples=false \
  -Dwith_gpulite=disabled -Duse_openmp=false -Deh_shared=true
meson compile -C "$lc_build"
meson install -C "$lc_build"

pc_dir=$prefix/lib/pkgconfig
if [ -n "${PKG_CONFIG_PATH-}" ]; then
  pc_dir=$pc_dir:$PKG_CONFIG_PATH
fi
export PKG_CONFIG_PATH=$pc_dir
pkg-config --print-errors --exists 'linkcell >= 0.6.0'
pkg-config --libs linkcell > "$base/linkcell-libs.txt"
grep -q -- '-llinkcell_eh' "$base/linkcell-libs.txt"

rm -rf "$seams_build"
meson setup "$seams_build" "$root" \
  --wrap-mode=nodownload --buildtype=debug \
  -Dwith_tests=true -Dwith_cli=false \
  -Dwith_chemfiles=disabled -Dwith_lammps=disabled \
  -Dwith_gromacs=disabled -Dwith_openmm=disabled \
  -Dwith_plumed=disabled -Dwith_gpulite=disabled \
  -Dwith_readcon=disabled -Dwith_minimage=disabled \
  -Dwith_blas=disabled -Dwith_ira=disabled \
  -Dwith_spglib=disabled -Dwith_nauty=disabled \
  -Dwith_sphericart=disabled -Dwith_backward=disabled \
  -Dwith_openmp=disabled

log=$seams_build/meson-logs/meson-log.txt
if grep -q "Executing subproject linkcell" "$log"; then
  echo "linkcell was configured from the subproject" >&2
  exit 1
fi
if ! grep -q "Determining dependency 'linkcell' with pkg-config" "$log"; then
  echo "linkcell was not resolved with pkg-config" >&2
  exit 1
fi
if ! grep -q "Run-time dependency linkcell found: YES" "$log"; then
  echo "installed linkcell was not found" >&2
  exit 1
fi

meson compile -C "$seams_build" exception_typeinfo

for bin in \
  "$seams_build/exception_typeinfo" \
  "$seams_build/libexception_typeinfo_hidden.so"
do
  readelf -d "$bin" > "$base/needed.txt"
  if ! grep -q 'liblinkcell_eh\.so' "$base/needed.txt"; then
    echo "$bin does not NEEDED liblinkcell_eh" >&2
    cat "$base/needed.txt" >&2
    exit 1
  fi
done

env -u LD_LIBRARY_PATH meson test -C "$seams_build" exception_typeinfo --print-errorlogs
echo "external linkcell_eh check ok"
