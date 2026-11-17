#!/bin/bash
# Copyright (c) d-SEAMS developers
# Install the Linux archive the way a user would, on a machine that
# does not have the build tree: check the checksum, run the command,
# compile a C program with pkg-config, and build the published dseams
# crate against that pkg-config file.
set -euo pipefail

if [[ -z ${VERSION:-} ]]; then
  printf '%s\n' "VERSION is required" >&2
  exit 1
fi

if command -v apt-get; then
  export DEBIAN_FRONTEND=noninteractive
  apt-get update
  apt-get install -y --no-install-recommends build-essential pkg-config curl ca-certificates
elif command -v dnf; then
  dnf install -y gcc pkgconf-pkg-config
else
  printf '%s\n' "no supported package manager" >&2
  exit 1
fi

compiler=$(command -v cc || command -v gcc)
printf '%s\n' "compiler: $compiler"
libc_ver=$(ldd --version)
printf '%s\n' "libc: ${libc_ver%%$'\n'*}"

test -f dist/SHA256SUMS
(cd dist && sha256sum -c SHA256SUMS)

archive=dist/libseams-${VERSION}-linux-x86_64.tar.gz
test -f "$archive"
stage=$(mktemp -d)
tar -C "$stage" -xzf "$archive"
prefix=$stage/libseams-${VERSION}
test -d "$prefix"
test -f "$prefix/LICENSE"
test -f "$prefix/licenses/eigen/COPYING.MPL2"
test -f "$prefix/licenses/eigen/COPYING.BSD"
test -f "$prefix/licenses/argum/LICENSE"
test -f "$prefix/licenses/backward-cpp/LICENSE"
test -f "$prefix/licenses/gcc-runtime/COPYING3"
test -f "$prefix/licenses/gcc-runtime/COPYING.RUNTIME"

count=0
for pc in "$prefix"/lib/pkgconfig/*.pc; do
  count=$((count + 1))
  line=$(grep '^prefix=' "$pc")
  if [[ $line != 'prefix=${pcfiledir}/../..' ]]; then
    printf '%s\n' "$(basename "$pc") prefix is $line" >&2
    exit 1
  fi
done
if [[ $count -ne 3 ]]; then
  printf '%s\n' "expected 3 pkg-config files, found $count" >&2
  exit 1
fi

moved=$(mktemp -d)
mv "$prefix" "$moved/"
prefix=$moved/libseams-${VERSION}

got=$(env -u LD_LIBRARY_PATH "$prefix/bin/seams" --version)
printf '%s\n' "$got"
if [[ $got != "seams $VERSION" ]]; then
  printf '%s\n' "version line is '$got'" >&2
  exit 1
fi

unset SEAMS_LIBDIR
unset SEAMS_LIB_DIR
unset SEAMS_SOURCE_DIR
unset CARGO_FEATURE_BUNDLED
export PKG_CONFIG_PATH=$prefix/lib/pkgconfig
pkg-config --exists --print-errors seams
if [[ $(pkg-config --modversion seams) != "$VERSION" ]]; then
  printf '%s\n' "pkg-config modversion is not $VERSION" >&2
  exit 1
fi
flags=$(pkg-config --cflags --libs seams)
printf '%s\n' "$flags"
resolved=$(realpath -- "$(pkg-config --variable=libdir seams)")
if [[ $resolved != "$(realpath -- "$prefix/lib")" ]]; then
  printf '%s\n' "pkg-config libdir resolved to $resolved" >&2
  exit 1
fi

cdir=$(mktemp -d)
cat > "$cdir/consumer.c" <<'EOF'
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
"$compiler" "$cdir/consumer.c" -o "$cdir/consumer" $(pkg-config --cflags --libs seams)
c_out=$(env LD_LIBRARY_PATH="$prefix/lib" "$cdir/consumer")
printf '%s\n' "$c_out"
if [[ $c_out != "c-consumer ok" ]]; then
  printf '%s\n' "C consumer failed" >&2
  exit 1
fi

export CARGO_HOME=$stage/cargo
export RUSTUP_HOME=$stage/rustup
curl --proto '=https' --tlsv1.2 -fsS https://sh.rustup.rs -o "$stage/rustup-init"
sh "$stage/rustup-init" -y --default-toolchain 1.85.0 --profile minimal
# rustup writes this file.
. "$CARGO_HOME/env"

rdir=$(mktemp -d)
mkdir -p "$rdir/src"
cat > "$rdir/Cargo.toml" <<'EOF'
[package]
name = "readme_example"
version = "0.1.0"
edition = "2021"

[dependencies]
dseams = "=2.14.0"
EOF
cat > "$rdir/src/main.rs" <<'EOF'
fn main() -> Result<(), seams::Error> {
    let a = 2.75 * 4.0 / 3.0_f64.sqrt();
    let sites = [
        (0.0, 0.0, 0.0), (0.5, 0.5, 0.0), (0.5, 0.0, 0.5), (0.0, 0.5, 0.5),
        (0.25, 0.25, 0.25), (0.75, 0.75, 0.25), (0.75, 0.25, 0.75), (0.25, 0.75, 0.75),
    ];
    let mut xyz = Vec::new();
    for (x, y, z) in sites {
        xyz.extend([x * a, y * a, z * a]);
    }
    let mut sys = seams::System::new()?;
    sys.set_ortho(&xyz, a, a, a)?;
    let ice = sys.chill_plus(&seams::ChillOptions::default())?;
    assert_eq!(ice.labels.len(), 8);
    Ok(())
}
EOF
(
  cd "$rdir"
  cargo build
)
example=$rdir/target/debug/readme_example
ldd_out=$(env LD_LIBRARY_PATH="$prefix/lib" ldd "$example")
printf '%s\n' "$ldd_out"
case $ldd_out in
  *"$prefix/lib/libseams.so"*) ;;
  *)
    printf '%s\n' "the example is not linked to the archive libseams" >&2
    exit 1
    ;;
esac
env LD_LIBRARY_PATH="$prefix/lib" "$example"
printf '%s\n' "dseams readme example ok"
