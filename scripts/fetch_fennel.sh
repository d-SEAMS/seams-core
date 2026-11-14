#!/bin/sh
# Fetch the Fennel compiler, one Lua file, for meson -Dfennel_lua.
# conda-forge does not package Fennel. Prints the file path.
# usage: scripts/fetch_fennel.sh DEST
set -eu

dest=${1:?usage: fetch_fennel.sh DEST}
version=1.5.3
sum=c6085d011692b8909ebb57088ad0f34da5f720b1d5ade170741335e1737371da
file="$dest/fennel-$version.lua"
mkdir -p "$dest"
if [ ! -f "$file" ] || ! printf '%s  %s\n' "$sum" "$file" | sha256sum -c - >/dev/null 2>&1; then
  curl -fsSL -o "$file.part" "https://fennel-lang.org/downloads/fennel-$version.lua"
  if ! printf '%s  %s\n' "$sum" "$file.part" | sha256sum -c - >/dev/null 2>&1; then
    rm -f "$file.part"
    echo "fennel-$version.lua does not match its sha256" >&2
    exit 1
  fi
  mv "$file.part" "$file"
fi
(cd "$dest" && printf '%s/fennel-%s.lua\n' "$(pwd)" "$version")
