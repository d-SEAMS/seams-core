#!/bin/sh
# Copy the engine sources that libseams_lang needs into src/seams, and the
# example dumps into inst/extdata, so that `R CMD build` makes a source
# package which compiles with no seams-core checkout and no network.
#
#   sh tools/vendor.sh [SEAMS_CORE_ROOT]
#
# The root defaults to $DSEAMS_ROOT, then to the checkout this package sits
# in (bindings/r/dseams). The list of engine objects is the seams/*.o block
# of src/Makevars.in, so that file is the one place to add a source.
set -e

pkg=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
root=${1:-${DSEAMS_ROOT:-$pkg/../../..}}
if [ ! -f "$root/src/mol_sys.cpp" ] || [ ! -f "$root/bindings/cxx/seams_lang.cpp" ]; then
  echo "vendor.sh: $root is not a seams-core checkout" >&2
  exit 1
fi
root=$(CDPATH= cd -- "$root" && pwd)

dest=$pkg/src/seams
rm -rf "$dest"
mkdir -p "$dest/include"

engine=$(sed -n 's|^[[:space:]]*seams/\([A-Za-z0-9_]*\)\.o.*|\1|p' "$pkg/src/Makevars.in")
if [ -z "$engine" ]; then
  echo "vendor.sh: no seams/*.o objects listed in src/Makevars.in" >&2
  exit 1
fi
engine_win=$(sed -n 's|^[[:space:]]*seams/\([A-Za-z0-9_]*\)\.o.*|\1|p' "$pkg/src/Makevars.win")
if [ "$engine" != "$engine_win" ]; then
  echo "vendor.sh: src/Makevars.in and src/Makevars.win list different engine objects" >&2
  exit 1
fi
version=$(sed -n 's/^Version:[[:space:]]*//p' "$pkg/DESCRIPTION")
if ! grep -q "SEAMS_LANG_VERSION='\"$version\"'" "$pkg/src/Makevars.win"; then
  echo "vendor.sh: src/Makevars.win does not carry version $version" >&2
  exit 1
fi
for name in $engine; do
  if [ "$name" = seams_lang ]; then
    cp "$root/bindings/cxx/seams_lang.cpp" "$dest/"
  else
    cp "$root/src/$name.cpp" "$dest/"
  fi
done
# seams_lang.cpp and seams_api.cpp share one shared object here, and both
# define seams_last_error and seams_path_nframes (with different
# signatures). The copied seams_lang.h renames the seams_lang pair, so
# dseams.cpp reaches seams_lang and seams_call.c reaches the C API.
awk '{ print } /^#define SEAMS_LANG_H_$/ {
  print "#define seams_last_error seams_lang_last_error"
  print "#define seams_path_nframes seams_lang_path_nframes"
}' "$root/bindings/cxx/seams_lang.h" > "$dest/seams_lang.h"
defs() {
  sed -n 's/^[A-Za-z_][A-Za-z0-9_ *]*[ *]\(seams_[A-Za-z0-9_]*\)(.*/\1/p' "$1" | sort -u
}
renamed=$(sed -n 's/^#define \(seams_[A-Za-z0-9_]*\) seams_lang_.*/\1/p' "$dest/seams_lang.h")
defs "$dest/seams_lang.cpp" | grep -vxF "$renamed" > "$dest/.lang_defs" || true
defs "$dest/seams_api.cpp" > "$dest/.api_defs"
clash=$(comm -12 "$dest/.lang_defs" "$dest/.api_defs")
rm -f "$dest/.lang_defs" "$dest/.api_defs"
if [ -n "$clash" ]; then
  echo "vendor.sh: seams_lang.cpp and seams_api.cpp both define:" $clash >&2
  exit 1
fi
cp "$root"/src/include/internal/*.h "$root"/src/include/internal/*.hpp "$dest/include/"
cp "$root/LICENSE" "$dest/LICENSE"

# Two edits make the copy an R package. [[nodiscard]] goes, as the meson
# build already silences it with -Wno-unused-result (R CMD check reports
# those warnings and a package may not pass -Wno-* flags). Console output
# and std::exit go through src/dseams_r_io.h (Writing R Extensions,
# "Printing" and "Writing portable packages").
for f in "$dest"/*.cpp "$dest"/include/*.h "$dest"/include/*.hpp; do
  sed -e 's/\[\[nodiscard\]\] *//g' "$f" > "$f.tmp"
  if grep -Eq 'std::(cerr|cout|exit)\b|fprintf\(stderr,' "$f.tmp"; then
    {
      echo '#include "dseams_r_io.h"'
      sed -e 's/std::cerr\b/dseams_r::err()/g' \
        -e 's/std::cout\b/dseams_r::out()/g' \
        -e 's/std::fprintf(stderr, */dseams_r::eprintf(/g' \
        -e 's/std::exit(/dseams_r::fail(/g' "$f.tmp"
    } > "$f"
    rm -f "$f.tmp"
  else
    mv "$f.tmp" "$f"
  fi
done
left=$(grep -En 'std::(cerr|cout|exit|abort|printf)\b|(^|[^a-z_:])(printf|puts)\(|\bstdout\b|fprintf\(stderr' \
  "$dest"/*.cpp "$dest"/include/*.h "$dest"/include/*.hpp || true)
if [ -n "$left" ]; then
  echo "vendor.sh: console calls left in the engine copy:" >&2
  echo "$left" >&2
  exit 1
fi

# The public C header for src/seams_call.c and seams_api.cpp. It goes in
# after the console check: its usage comment prints, and nothing in the
# package compiles that comment.
cp "$root/include/seams.h" "$dest/include/"
major=${version%%.*}
rest=${version#*.}
minor=${rest%%.*}
patch=${rest#*.}
sed -e "s/@MAJOR@/$major/" -e "s/@MINOR@/$minor/" -e "s/@PATCH@/$patch/" \
  -e "s/@STRING@/$version/" "$root/include/seams_version.h.in" > "$dest/include/seams_version.h"

rev=unknown
if command -v git >/dev/null 2>&1 && git -C "$root" rev-parse HEAD >/dev/null 2>&1; then
  rev=$(git -C "$root" rev-parse HEAD)
  if ! git -C "$root" diff --quiet HEAD -- src include bindings/cxx; then
    rev="$rev-dirty"
  fi
fi
{
  echo "d-SEAMS engine sources, copied by tools/vendor.sh."
  echo "seams-core commit: $rev"
  echo "Upstream: https://github.com/d-SEAMS/seams-core"
} > "$dest/README"

extdata=$pkg/inst/extdata
mkdir -p "$extdata"
for f in ice1h.lammpstrj ice1h.source.txt exampleTraj.lammpstrj \
  nucleation_f1.lammpstrj nucleation_f1.source.txt; do
  cp "$root/input/traj/$f" "$extdata/"
done

echo "vendor.sh: $(echo $engine | wc -w | tr -d ' ') engine sources from $rev"
