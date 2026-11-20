#!/bin/sh
# Copyright (c) d-SEAMS developers
# Unpack chemfiles' molfiles.tar.gz, apply the reader fixes, and repack it.
# chemfiles' external/CMakeLists.txt unpacks the archive when it configures,
# and it refuses to unpack again while the archive's sha256 is unchanged.
# A second run of this script sees the patch already applied and leaves the
# archive untouched, so that hash stays put.
set -eu
here=$(CDPATH= cd -- "$(dirname "$0")" && pwd)
patch_file=$here/chemfiles-0.10.4-molfile-sources.patch
if [ ! -f external/molfiles.tar.gz ]; then
  echo "chemfiles molfiles archive is missing" >&2
  exit 1
fi
if [ ! -f "$patch_file" ]; then
  echo "molfile source patch is missing: $patch_file" >&2
  exit 1
fi
tmp=$(mktemp -d)
trap 'rm -rf "$tmp"' EXIT
tar -xzf external/molfiles.tar.gz -C "$tmp"
if patch -p1 -R -s -f --dry-run -d "$tmp" -i "$patch_file" >/dev/null 2>&1; then
  exit 0
fi
patch -p1 -d "$tmp" -i "$patch_file"
tar -cf - -C "$tmp" molfiles | gzip -n > external/molfiles.tar.gz
