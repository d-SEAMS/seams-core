#!/bin/sh
# Copyright (c) d-SEAMS developers
# chemfiles 0.10.4 vendors zlib 1.2.11. That zutil.h treats
# TARGET_OS_MAC as classic Mac OS and replaces fdopen with NULL.
# Apple's SDK defines TARGET_OS_MAC on macOS, which does have fdopen,
# and the stub then breaks the stdio.h declaration. Classic Mac OS is
# the MACOS macro alone. __APPLE__ still sets OS_CODE later in the header.
set -eu
if [ ! -f external/zlib.tar.gz ]; then
  exit 0
fi
tmp=$(mktemp -d)
trap 'rm -rf "$tmp"' EXIT
tar -xzf external/zlib.tar.gz -C "$tmp"
header=$tmp/zlib/zutil.h
if [ ! -f "$header" ]; then
  echo "chemfiles zlib archive has no zutil.h" >&2
  exit 1
fi
if grep -q 'defined(MACOS) || defined(TARGET_OS_MAC)' "$header"; then
  sed 's/defined(MACOS) || defined(TARGET_OS_MAC)/defined(MACOS)/' "$header" > "$header.new"
  mv "$header.new" "$header"
  tar -czf external/zlib.tar.gz -C "$tmp" zlib
fi
