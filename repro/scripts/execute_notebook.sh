#!/usr/bin/env bash
# Convert a percent-format notebook and execute it with papermill.
set -euo pipefail
if [[ $# -ne 2 ]]; then
  echo "usage: $0 <notebook.py> <out.ipynb>" >&2
  exit 2
fi
src=$1
out=$2
mkdir -p "$(dirname "$out")"
tmp=$(mktemp --suffix=.ipynb)
trap 'rm -f "$tmp"' EXIT
jupytext --to ipynb --output "$tmp" "$src"
papermill "$tmp" "$out" --cwd "$(pwd)" -k python3
