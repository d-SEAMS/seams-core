#!/bin/sh
# Copyright (c) d-SEAMS developers
# Check every patch under cmake/patches against pristine chemfiles 0.10.4.
#
# Usage: scripts/check_chemfiles_patches.sh [chemfiles-0.10.4.tar.gz]
# Without an argument the release tarball is downloaded. Its sha256 is
# checked either way.
#
# The patches go on in the order cmake/ChemfilesDep.cmake applies them.
# Each one must pass both `git apply --check` and `patch --dry-run` on the
# tree the earlier patches produced, and is then applied with patch. The
# *-sources patches are checked inside the external archive they edit.
# A patch file in cmake/patches that this list does not name is an error,
# so a new patch cannot go unchecked.
set -eu
root=$(CDPATH= cd -- "$(dirname "$0")/.." && pwd)
patches=$root/cmake/patches
sha256=b8232ddaae2953538274982838aa6c2df87d300f7e2f80e92c171581e06325ba
url=https://github.com/chemfiles/chemfiles/archive/refs/tags/0.10.4.tar.gz

tree_patches="mmcif mol2 frame-visit lammpsdata-section molfiles"
archive_patches="molfiles:molfile-sources tng:tng-sources mmtf-cpp:mmtf-sources gemmi:gemmi-sources"
late_patches="hostile span-unused"

work=$(mktemp -d)
trap 'rm -rf "$work"' EXIT

if [ $# -ge 1 ]; then
  tarball=$1
else
  tarball=$work/chemfiles-0.10.4.tar.gz
  curl -fsSL -o "$tarball" "$url"
fi
if command -v sha256sum >/dev/null 2>&1; then
  got=$(sha256sum "$tarball" | cut -d' ' -f1)
else
  got=$(shasum -a 256 "$tarball" | cut -d' ' -f1)
fi
if [ "$got" != "$sha256" ]; then
  echo "unexpected chemfiles tarball sha256 $got (want $sha256)" >&2
  exit 1
fi

listed=" $tree_patches $late_patches "
for pair in $archive_patches; do
  listed="$listed${pair#*:} "
done
status=0
for file in "$patches"/*.patch; do
  name=${file##*/chemfiles-0.10.4-}
  name=${name%.patch}
  case "$listed" in
    *" $name "*) ;;
    *) echo "FAIL $file is not in the checked list" >&2; status=1 ;;
  esac
done

# check NAME DIR: both tools must accept the patch on DIR, then apply it
check() {
  file=$patches/chemfiles-0.10.4-$1.patch
  if ! (cd "$2" && git apply --check -p1 "$file"); then
    echo "FAIL $1: git apply --check" >&2
    return 1
  fi
  if ! patch -p1 -N -s --dry-run -d "$2" -i "$file" >/dev/null; then
    echo "FAIL $1: patch --dry-run" >&2
    return 1
  fi
  patch -p1 -N -s -d "$2" -i "$file" >/dev/null
  echo "ok   $1"
}

tar -xzf "$tarball" -C "$work"
src=$work/chemfiles-0.10.4
for name in $tree_patches; do
  check "$name" "$src" || status=1
done
for pair in $archive_patches; do
  dir=${pair%%:*}
  unpacked=$work/archive-$dir
  mkdir -p "$unpacked"
  tar -xzf "$src/external/$dir.tar.gz" -C "$unpacked"
  check "${pair#*:}" "$unpacked" || status=1
done
for name in $late_patches; do
  check "$name" "$src" || status=1
done
# The build's patch step ends with these two archive edits.
for script in zlib-apple fmt-libcxx; do
  if (cd "$src" && sh "$patches/chemfiles-0.10.4-$script.sh"); then
    echo "ok   $script.sh"
  else
    echo "FAIL $script.sh" >&2
    status=1
  fi
done
exit $status
