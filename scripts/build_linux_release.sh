#!/bin/bash
# Copyright (c) d-SEAMS developers
# Build dist/libseams-<version>-linux-x86_64.tar.gz from the source
# archive, inside a manylinux_2_28 image. The image's libc is glibc
# 2.28 and its compiler leaves newer C++ runtime objects inside the
# binaries, so the archive runs where glibc is 2.28 or newer.
set -euo pipefail

root=$(CDPATH= cd -- "$(dirname "$0")/.." && pwd)
cd "$root"

version=$(python3 ci/version/seams_version.py)
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
# Root would keep the uid stored in the archive, and git would refuse
# the readcon checkout.
tar --no-same-owner -C "$work/src" -xzf "$src_tar"
src=$work/src/seams-core-${version}
# readcon-core v0.16.2, readcon-db v0.1.8. nodownload will not clone them.
# The dist tree already has the wrap checkout, with no .git. Replace it.
pin_repo() {
  local url=$1 sha=$2 dest=$3 got
  rm -rf "$dest"
  git init "$dest"
  git -C "$dest" remote add origin "$url"
  git -C "$dest" fetch --depth 1 origin "$sha"
  git -C "$dest" checkout --detach FETCH_HEAD
  got=$(git -C "$dest" rev-parse HEAD)
  if [[ $got != "$sha" ]]; then
    printf '%s\n' "checkout $got is not $sha" >&2
    exit 1
  fi
}
pin_repo https://github.com/lode-org/readcon-core.git \
  9e12ea9727d5f3c85d8eeb4d927beb634e1bdb59 \
  "$src/subprojects/readcon-core"
pin_repo https://github.com/lode-org/readcon-db.git \
  0001a2c5a931bdde6bc8fe9b308f666dc28f4fe5 \
  "$src/subprojects/readcon-db"
export CARGO_HOME=$work/cargo
export RUSTUP_HOME=$work/rustup
# rustup 1.29.3
rustup_sha256=b1c8cb14e30889e7e2999e6b69e62d9e9b9ffc0f4df8180564d306d9680e6c55
curl --proto '=https' --tlsv1.2 -fsS \
  -o "$work/rustup-init" \
  https://static.rust-lang.org/rustup/archive/1.29.3/x86_64-unknown-linux-gnu/rustup-init
echo "${rustup_sha256}  $work/rustup-init" | sha256sum -c -
chmod +x "$work/rustup-init"
"$work/rustup-init" -y --default-toolchain 1.88.0 --profile minimal
# rustup writes this file.
. "$CARGO_HOME/env"
printf '%s\n' "rustc: $(rustc --version)"
prefix=$work/libseams-${version}
export PREFIX=$prefix
export CPU_COUNT=${CPU_COUNT:-$(nproc)}
export SEAMS_WITH_READCON=enabled
(cd "$src" && bash recipe/build-libseams.sh)

EIGEN_SRC=$eigen_src \
  READCON_CORE_SRC=$src/subprojects/readcon-core \
  READCON_DB_SRC=$src/subprojects/readcon-db \
  bash "$root/scripts/finish_libseams_prefix.sh" "$prefix"

export PKG_CONFIG_PATH=$prefix/lib/pkgconfig
pkg-config --exists --print-errors seams
resolved=$(realpath -- "$(pkg-config --variable=libdir seams)")
if [[ $resolved != "$(realpath -- "$prefix/lib")" ]]; then
  printf '%s\n' "pkg-config libdir resolved to $resolved" >&2
  exit 1
fi
printf '%s\n' "pkg-config: $(pkg-config --cflags --libs seams)"

python3 - "$prefix" <<'PY' | tee "$work/symbol-scan.txt"
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
needs_gomp = False
elfs = [path for path in prefix.joinpath("bin").iterdir() if path.is_file()]
elfs += list(prefix.joinpath("lib").glob("*.so"))
for elf in sorted(elfs):
    dynamic = readelf(elf, "-d")
    for line in dynamic.splitlines():
        if "NEEDED" in line and "[libgomp.so.1]" in line:
            needs_gomp = True
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
print("libgomp=" + ("yes" if needs_gomp else "no"))
PY
glibc=$(sed -n 's/^GLIBC_//p' "$work/symbol-scan.txt")
gomp=$(sed -n 's/^libgomp=//p' "$work/symbol-scan.txt")
if [[ -z $glibc || -z $gomp ]]; then
  printf '%s\n' "symbol scan did not report glibc and libgomp" >&2
  exit 1
fi
dseams_version=$(awk -F= '$1 == "dseams" { print $2 }' "$root/scripts/linux_tarball_consumers")
dseams_fortran=$(awk -F= '$1 == "dseams_fortran" { print $2 }' "$root/scripts/linux_tarball_consumers")
if [[ -z $dseams_version || -z $dseams_fortran ]]; then
  printf '%s\n' "binding version file is incomplete" >&2
  exit 1
fi
pc_line=$(grep '^prefix=' "$prefix/lib/pkgconfig/seams.pc")
{
  printf '%s\n' "The linux-x86_64 archive from this build needs glibc ${glibc} or newer."
  if [[ $gomp == yes ]]; then
    printf '%s\n' "It loads libgomp.so.1 from the system. Debian and Ubuntu ship that runtime as libgomp1, and RHEL ships it as libgomp."
  fi
  printf '%s\n' "seams.pc records ${pc_line}."
  if [[ -f $prefix/LICENSE && -d $prefix/licenses ]]; then
    printf '%s\n' "The archive includes LICENSE and the licence texts of the code compiled into the binaries."
  fi
  if [[ -f $prefix/lib/libreadcon_core.so ]]; then
    printf '%s\n' "The .con reader is included."
  fi
  printf '%s\n' "The release test installs dseams ${dseams_version} from crates.io and dseams-fortran ${dseams_fortran}."
} > "$out/linux-archive-note.md"
printf '%s\n' "linux archive note:"
cat "$out/linux-archive-note.md"

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

cat > "$work/con_reader.c" <<'EOF'
#include "seams.h"
#include <stdio.h>
int main(int argc, char **argv) {
  seams_system *sys = seams_system_create();
  if (argc != 2 || sys == 0) {
    return 1;
  }
  if (seams_system_read(sys, argv[1], 1, -1) != SEAMS_OK) {
    return 1;
  }
  if (seams_system_count(sys) <= 0) {
    return 1;
  }
  seams_system_destroy(sys);
  printf("con-reader ok\n");
  return 0;
}
EOF
gcc "$work/con_reader.c" -o "$work/con_reader" $(pkg-config --cflags --libs seams)
env LD_LIBRARY_PATH="$prefix/lib" "$work/con_reader" \
  "$src/subprojects/readcon-core/resources/test/tiny_multi_cuh2.con"

mkdir -p "$out"
tar -C "$work" -czf "$out/libseams-${version}-linux-x86_64.tar.gz" "libseams-${version}"
printf '%s\n' "wrote libseams-${version}-linux-x86_64.tar.gz"
