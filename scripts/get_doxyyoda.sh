#!/usr/bin/env bash
# Fetch the doxyYoda HTML theme next to Doxyfile-prj.cfg.
set -euo pipefail
ver=0.2.2
root=$(git rev-parse --show-toplevel)
cd "$root"
if [ -d doxyYoda/html ]; then
  echo "doxyYoda already present"
  exit 0
fi
curl -fSL "https://github.com/HaoZeke/doxyYoda/releases/download/v${ver}/doxyYoda_${ver}.tar.gz" -o doxyYoda.tar.gz
tar xf doxyYoda.tar.gz
rm -f doxyYoda.tar.gz
echo "extracted doxyYoda ${ver}"
