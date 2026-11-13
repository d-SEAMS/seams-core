#!/bin/sh
# Re-run the suite under AddressSanitizer and UndefinedBehaviorSanitizer.
#
#   scripts/sanitize_suite.sh
#   scripts/sanitize_suite.sh --reconfigure
#
# Extra arguments are passed to `meson setup`. The build directory is
# build-san next to this repository. ThreadSanitizer is a separate
# configure: -Db_sanitize=thread.
set -eu
root=$(CDPATH= cd -- "$(dirname "$0")/.." && pwd)
build=$root/build-san
reconfigure=
if [ "${1:-}" = "--reconfigure" ]; then
  reconfigure=--reconfigure
  shift
fi
if [ -d "$build" ] && [ -z "$reconfigure" ]; then
  reconfigure=--reconfigure
fi
meson setup "$build" $reconfigure \
  --buildtype=debugoptimized \
  -Db_sanitize=address,undefined \
  -Dwith_tests=true \
  -Dwith_python=true \
  -Dwith_cli=true \
  "$@"
ninja -C "$build"
meson test -C "$build" --print-errorlogs
