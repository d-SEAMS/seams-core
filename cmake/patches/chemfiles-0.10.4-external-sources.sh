#!/bin/sh
# Copyright (c) d-SEAMS developers
# Usage: chemfiles-0.10.4-external-sources.sh DIR NAME
# Unpack chemfiles' external/DIR.tar.gz, apply chemfiles-0.10.4-NAME.patch
# next to this script, and repack the archive.
# chemfiles' external/CMakeLists.txt unpacks the archive when it configures,
# and it refuses to unpack again while the archive's sha256 is unchanged.
# A second run of this script sees the patch already applied and leaves the
# archive untouched, so that hash stays put.
set -eu
if [ $# -ne 2 ]; then
  echo "usage: $0 DIR NAME" >&2
  exit 2
fi
dir=$1
here=$(CDPATH= cd -- "$(dirname "$0")" && pwd)
patch_file=$here/chemfiles-0.10.4-$2.patch
archive=external/$dir.tar.gz
if [ ! -f "$archive" ]; then
  echo "chemfiles archive is missing: $archive" >&2
  exit 1
fi
if [ ! -f "$patch_file" ]; then
  echo "source patch is missing: $patch_file" >&2
  exit 1
fi
tmp=$(mktemp -d)
trap 'rm -rf "$tmp"' EXIT
tar -xzf "$archive" -C "$tmp"
if patch -p1 -R -s -f --dry-run -d "$tmp" -i "$patch_file" >/dev/null 2>&1; then
  exit 0
fi
patch -p1 -d "$tmp" -i "$patch_file"
tar -cf - -C "$tmp" "$dir" | gzip -n > "$archive"
