#!/bin/sh
# Fetch the class headers of the LAMMPS release that lmp (or $LMP) reports,
# for meson -Dlammps_src_include. Prints the src directory.
# usage: scripts/fetch_lammps_headers.sh DEST
set -eu

dest=${1:?usage: fetch_lammps_headers.sh DEST}
lmp=${LMP:-lmp}
banner=$(echo "info config" | "$lmp" -log none -screen none -echo none 2>/dev/null; "$lmp" -h 2>/dev/null | head -n 5)
version=$(printf '%s\n' "$banner" | sed -n 's/.*- \([0-9][0-9]* [A-Z][a-z][a-z] [0-9]\{4\}\)\( - Update \([0-9]*\)\)\{0,1\}$/\1|\3/p' | head -n 1)
if [ -z "$version" ]; then
  version=$(printf '%s\n' "$banner" | sed -n 's/^LAMMPS (\([0-9][0-9]* [A-Z][a-z][a-z] [0-9]\{4\}\)\( - Update \([0-9]*\)\)\{0,1\}).*/\1|\3/p' | head -n 1)
fi
if [ -z "$version" ]; then
  echo "cannot read the LAMMPS version from $lmp" >&2
  exit 1
fi
date=${version%|*}
update=${version#*|}
set -- $date
day=${1#0} month=$2 year=$3
stamp="${day}${month}${year}"
if [ -n "$update" ]; then
  tags="stable_${stamp}_update${update}"
else
  tags="stable_${stamp} patch_${stamp}"
fi
for tag in $tags; do
  rm -rf "$dest"
  if git clone -q --depth 1 --branch "$tag" --filter=blob:none --sparse \
      https://github.com/lammps/lammps.git "$dest" 2>/dev/null; then
    git -C "$dest" sparse-checkout set --no-cone '/src/*.h' '/src/fmt/' '/src/STUBS/' >/dev/null
    (cd "$dest/src" && pwd)
    exit 0
  fi
done
echo "no LAMMPS tag among: $tags" >&2
exit 1
