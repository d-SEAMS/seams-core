#!/usr/bin/env bash
# Fetch the doxyYoda HTML theme into docs/doxyYoda/.
set -euo pipefail
ver=0.2.2
script_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
root=$(cd -- "$script_dir/.." && pwd)
dest="$root/docs/doxyYoda"
if [ -f "$dest/html/header.html" ]; then
  echo "doxyYoda already present"
  exit 0
fi
tmp=$(mktemp)
trap 'rm -f "$tmp"' EXIT
curl -fSL "https://github.com/HaoZeke/doxyYoda/releases/download/v${ver}/doxyYoda_${ver}.tar.gz" -o "$tmp"
mkdir -p "$root/docs"
tar -xf "$tmp" -C "$root/docs"
echo "extracted doxyYoda ${ver} -> docs/doxyYoda"
