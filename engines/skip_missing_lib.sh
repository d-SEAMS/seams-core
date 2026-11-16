#!/bin/sh
# Copyright (c) d-SEAMS developers
# Run an engine program. Exit 77 when its shared library is not installed.
set -u
if [ "$#" -lt 1 ]; then
  echo "usage: skip_missing_lib.sh PROGRAM [ARGS...]" >&2
  exit 2
fi
# Under the sanitizer suite LD_PRELOAD holds libasan. mktemp, cat, grep and
# rm are not seams binaries and fail LeakSanitizer at exit, so only the
# engine program gets the preload.
preload=${LD_PRELOAD-}
unset LD_PRELOAD
err=$(mktemp)
out=$(mktemp)
if [ -n "$preload" ]; then
  LD_PRELOAD=$preload "$@" >"$out" 2>"$err"
else
  "$@" >"$out" 2>"$err"
fi
status=$?
if [ "$status" -eq 0 ]; then
  cat "$out"
  rm -f "$out" "$err"
  exit 0
fi
cat "$err" >&2
if grep -q "cannot open shared object file" "$err" || [ "$status" -eq 127 ]; then
  echo "engine skipped: shared library was not found" >&2
  rm -f "$out" "$err"
  exit 77
fi
rm -f "$out" "$err"
exit "$status"
