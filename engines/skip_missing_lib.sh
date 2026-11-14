#!/bin/sh
# Run an engine program. Exit 77 when its shared library is not installed.
set -u
if [ "$#" -lt 1 ]; then
  echo "usage: skip_missing_lib.sh PROGRAM [ARGS...]" >&2
  exit 2
fi
err=$(mktemp)
out=$(mktemp)
if "$@" >"$out" 2>"$err"; then
  cat "$out"
  rm -f "$out" "$err"
  exit 0
fi
status=$?
cat "$err" >&2
if grep -q "cannot open shared object file" "$err" || [ "$status" -eq 127 ]; then
  echo "engine skipped: shared library was not found" >&2
  rm -f "$out" "$err"
  exit 77
fi
rm -f "$out" "$err"
exit "$status"
