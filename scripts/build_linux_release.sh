#!/bin/bash
# Copyright (c) d-SEAMS developers
# Build dist/libseams-<version>-linux-x86_64.tar.gz from the source
# archive, inside a manylinux_2_28 image. The image's libc is glibc
# 2.28 and its compiler leaves newer C++ runtime objects inside the
# binaries, so the archive runs where glibc is 2.28 or newer.
set -euo pipefail

root=$(CDPATH= cd -- "$(dirname "$0")/.." && pwd)
cd "$root"

version=$(python3 scripts/check_release_version.py --print)
if [[ -n ${VERSION:-} && $VERSION != "$version" ]]; then
  printf '%s\n' "VERSION=$VERSION does not match meson.build ($version)" >&2
  exit 1
fi

src_tar=$root/dist/seams-core-${version}.tar.gz
if [[ ! -f $src_tar ]]; then
  printf '%s\n' "missing $src_tar" >&2
  exit 1
fi
out=${OUT_DIR:-$root/dist}
mkdir -p "$out"

cc_ver=$(g++ --version)
libc_ver=$(ldd --version)
printf '%s\n' "compiler: ${cc_ver%%$'\n'*}"
printf '%s\n' "libc: ${libc_ver%%$'\n'*}"

work=$(mktemp -d)
cleanup() { rm -rf "$work"; }
trap cleanup EXIT

python3 -m venv "$work/venv"
"$work/venv/bin/pip" install --disable-pip-version-check --upgrade pip
"$work/venv/bin/pip" install --disable-pip-version-check meson ninja
export PATH="$work/venv/bin:$PATH"
printf '%s\n' "meson: $(meson --version)"

python3 - "$work" > "$work/eigen-info.txt" <<'PY'
import hashlib
import pathlib
import sys
import tarfile
import urllib.request

work = pathlib.Path(sys.argv[1])
vals = {}
for line in pathlib.Path("subprojects/eigen.wrap").read_text().splitlines():
    if "=" not in line or line.startswith("["):
        continue
    key, value = line.split("=", 1)
    vals[key.strip()] = value.strip()
urls = [vals["source_url"]]
if vals.get("source_fallback_url"):
    urls.append(vals["source_fallback_url"])
blob = work / "eigen.tar.bz2"
digest = vals["source_hash"]
last = "no eigen url"
for url in urls:
    try:
        urllib.request.urlretrieve(url, blob)
    except Exception as exc:
        last = str(exc)
        continue
    got = hashlib.sha256(blob.read_bytes()).hexdigest()
    if got == digest:
        break
    last = "hash " + got
else:
    raise SystemExit("eigen download failed: " + last)
unpack = work / "eigen-unpack"
unpack.mkdir()
with tarfile.open(blob, "r:bz2") as archive:
    archive.extractall(unpack, filter="data")
src = unpack / vals["directory"]
macros = (src / "Eigen/src/Core/util/Macros.h").read_text()

def field(name):
    for line in macros.splitlines():
        parts = line.split()
        if len(parts) >= 3 and parts[0] == "#define" and parts[1] == name:
            return parts[2]
    raise SystemExit("missing " + name)

version = ".".join(
    field(name)
    for name in ("EIGEN_WORLD_VERSION", "EIGEN_MAJOR_VERSION", "EIGEN_MINOR_VERSION")
)
print(src)
print(version)
PY
eigen_src=$(sed -n '1p' "$work/eigen-info.txt")
eigen_version=$(sed -n '2p' "$work/eigen-info.txt")
printf '%s\n' "eigen: $eigen_version"

eigen_prefix=$work/eigen-prefix
mkdir -p "$eigen_prefix/include/eigen3" "$eigen_prefix/share/pkgconfig"
cp -a "$eigen_src/Eigen" "$eigen_prefix/include/eigen3/"
if [[ -d $eigen_src/unsupported ]]; then
  cp -a "$eigen_src/unsupported" "$eigen_prefix/include/eigen3/"
fi
if [[ -f $eigen_src/signature_of_eigen3_matrix_library ]]; then
  cp -a "$eigen_src/signature_of_eigen3_matrix_library" "$eigen_prefix/include/eigen3/"
fi
cat > "$eigen_prefix/share/pkgconfig/eigen3.pc" <<EOF
prefix=${eigen_prefix}
includedir=\${prefix}/include/eigen3
Name: Eigen3
Description: A C++ template library for linear algebra
Version: ${eigen_version}
Cflags: -I\${includedir}
EOF
export PKG_CONFIG_PATH=$eigen_prefix/share/pkgconfig
pkg-config --exists --print-errors "eigen3 >= 3.4.0"

mkdir -p "$work/src"
tar -C "$work/src" -xzf "$src_tar"
src=$work/src/seams-core-${version}
prefix=$work/libseams-${version}
export PREFIX=$prefix
export CPU_COUNT=${CPU_COUNT:-$(nproc)}
(cd "$src" && bash recipe/build-libseams.sh)

EIGEN_SRC=$eigen_src bash "$root/scripts/finish_libseams_prefix.sh" "$prefix"

export PKG_CONFIG_PATH=$prefix/lib/pkgconfig
pkg-config --exists --print-errors seams
resolved=$(realpath -- "$(pkg-config --variable=libdir seams)")
if [[ $resolved != "$(realpath -- "$prefix/lib")" ]]; then
  printf '%s\n' "pkg-config libdir resolved to $resolved" >&2
  exit 1
fi
printf '%s\n' "pkg-config: $(pkg-config --cflags --libs seams)"

python3 - "$prefix" <<'PY'
import re
import subprocess
import sys
from pathlib import Path

prefix = Path(sys.argv[1])
name_re = re.compile(r"Name: ((?:GLIBC|GLIBCXX|CXXABI|GOMP|OMP)_[0-9.]+)")

def readelf(path, flag):
    done = subprocess.run(
        ["readelf", flag, str(path)],
        check=True,
        capture_output=True,
        text=True,
    )
    return done.stdout

def defined(path):
    text = readelf(path, "-V")
    start = text.find("Version definition section")
    if start < 0:
        return set()
    stop = text.find("Version needs section", start)
    body = text[start:] if stop < 0 else text[start:stop]
    return set(name_re.findall(body))

def needed(path):
    text = readelf(path, "-V")
    start = text.find("Version needs section")
    if start < 0:
        return set()
    return set(name_re.findall(text[start:]))

def parts(symbol):
    kind, ver = symbol.split("_", 1)
    return kind, tuple(int(piece) for piece in ver.split("."))

provided = set()
for soname in (
    "libc.so.6",
    "libm.so.6",
    "libdl.so.2",
    "libpthread.so.0",
    "libgcc_s.so.1",
    "libstdc++.so.6",
    "libgomp.so.1",
):
    located = subprocess.check_output(
        ["g++", "-print-file-name=" + soname], text=True
    ).strip()
    path = Path(located)
    if not path.is_file():
        raise SystemExit("compiler did not locate " + soname)
    provided |= defined(path)

highest = {}
elfs = [path for path in prefix.joinpath("bin").iterdir() if path.is_file()]
elfs += list(prefix.joinpath("lib").glob("*.so"))
for elf in sorted(elfs):
    dynamic = readelf(elf, "-d")
    for line in dynamic.splitlines():
        if "RPATH" in line or "RUNPATH" in line:
            found = re.search(r"\[(.+)\]", line)
            if found is None or found.group(1) not in {"$ORIGIN", "$ORIGIN/../lib"}:
                raise SystemExit("unexpected runpath in " + elf.name + ": " + line.strip())
    for symbol in sorted(needed(elf)):
        if symbol not in provided:
            raise SystemExit(elf.name + " needs " + symbol + ", missing from the image libc")
        kind, ver = parts(symbol)
        if kind == "GLIBC" and ver > (2, 28):
            raise SystemExit(elf.name + " needs GLIBC_" + ".".join(map(str, ver)))
        highest[kind] = max(highest.get(kind, ()), ver)

for kind, ver in sorted(highest.items()):
    print(kind + "_" + ".".join(str(piece) for piece in ver))
PY

got=$(env -u LD_LIBRARY_PATH "$prefix/bin/seams" --version)
printf '%s\n' "seams: $got"
if [[ $got != "seams $version" ]]; then
  printf '%s\n' "version line is '$got'" >&2
  exit 1
fi

cat > "$work/consumer.c" <<'EOF'
#include "seams.h"
#include <stdio.h>
int main(void) {
  seams_system *sys = seams_system_create();
  if (sys == 0) {
    return 1;
  }
  seams_system_destroy(sys);
  printf("c-consumer ok\n");
  return 0;
}
EOF
gcc "$work/consumer.c" -o "$work/consumer" $(pkg-config --cflags --libs seams)
env LD_LIBRARY_PATH="$prefix/lib" "$work/consumer"

mkdir -p "$out"
tar -C "$work" -czf "$out/libseams-${version}-linux-x86_64.tar.gz" "libseams-${version}"
printf '%s\n' "wrote libseams-${version}-linux-x86_64.tar.gz"
