#!/usr/bin/env bash
# Copy campaign outputs into repro/zenodo/payload/. No manuscript.
set -euo pipefail
ROOT=$(cd -- "$(dirname -- "$0")/.." && pwd)
SRC=$ROOT/repro/results
DST=$ROOT/repro/zenodo/payload
test -f "$SRC/paper_manifest.json" || {
  echo "missing $SRC/paper_manifest.json; run the campaign first" >&2
  exit 1
}
if [ -e "$DST" ]; then
  BACKUP=$DST.previous.$(date -u +%Y%m%dT%H%M%SZ)
  test ! -e "$BACKUP" || {
    echo "backup path exists: $BACKUP" >&2
    exit 1
  }
  mv "$DST" "$BACKUP"
  echo "preserved existing payload at $BACKUP"
fi
mkdir -p "$DST"
cp -a "$SRC/paper_manifest.json" "$DST/"
cp -a "$SRC/source-manifest.json" "$DST/"
cp -a "$SRC/workflow-parity.json" "$DST/"
cp -a "$SRC/conditions.txt" "$DST/"
for f in "$SRC"/tip-*.txt "$SRC"/base-*.txt "$SRC"/trajectory-incremental.txt \
         "$SRC"/ql-*.txt "$SRC"/figshare-incremental.json; do
  [ -f "$f" ] && cp -a "$f" "$DST/"
done
cp -a "$ROOT/repro/config.yaml" "$DST/"
cp -a "$ROOT/repro/ecosystem-lock.json" "$DST/"
cp -a "$ROOT/repro/Snakefile" "$DST/"
if [ -d "$SRC/figshare-demos" ]; then
  mkdir -p "$DST/figshare-demos"
  cp -a "$SRC/figshare-demos/"*.json "$DST/figshare-demos/" 2>/dev/null || true
fi
# Exclusive-node identity for GPU and v1 bars the paper cites
for f in \
  "$ROOT/repro/reference/gpu-conditions-elja-1787449.txt" \
  "$ROOT/repro/reference/tip-gpu-batch-elja-1787449.txt" \
  "$ROOT/repro/reference/conditions-elja-1785712.txt"
do
  [ -f "$f" ] && cp -a "$f" "$DST/"
done
(cd "$DST" && find . -type f ! -name MANIFEST | sed 's|^\./||' | sort > MANIFEST)
echo "staged $DST ($(wc -l < "$DST/MANIFEST") files)"
