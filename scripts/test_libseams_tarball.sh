#!/bin/bash
# Copyright (c) d-SEAMS developers
# Install the Linux archive the way a user would, on a machine that
# does not have the build tree. CHECK selects one of three runs:
#   c        compile docs/examples/quickstart.c with pkg-config,
#            then again with the archive rpath and no library path
#   rust     cargo add dseams@2.14.0, without the bundled feature
#   fortran  dseams-fortran v2.14.0 `fpm test`
set -euo pipefail

if [[ -z ${VERSION:-} || -z ${CHECK:-} ]]; then
  printf '%s\n' "VERSION and CHECK are required" >&2
  exit 1
fi

case $CHECK in
  c|rust)
    if command -v apt-get; then
      export DEBIAN_FRONTEND=noninteractive
      apt-get update
      apt-get install -y --no-install-recommends build-essential pkg-config ca-certificates
    else
      printf '%s\n' "no supported package manager" >&2
      exit 1
    fi
    ;;
  fortran)
    export DEBIAN_FRONTEND=noninteractive
    apt-get update
    apt-get install -y --no-install-recommends gfortran pkg-config ca-certificates curl git
    ;;
  *)
    printf '%s\n' "unknown CHECK=$CHECK" >&2
    exit 1
    ;;
esac

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
test -f "$prefix/licenses/readcon-core/LICENSE"
test -f "$prefix/licenses/readcon-db/LICENSE"
test -f "$prefix/lib/libreadcon_core.so"
test -f "$prefix/lib/libreadcon_db.so"
test -f "$prefix/lib/libyodaLib.so"
test -f "$prefix/lib/libseams.so"
test -f "$prefix/lib/libseams_lang.so"
test -f "$prefix/include/seams.h"
test -f "$prefix/include/seams_version.h"
test -f "$prefix/include/seams.hpp"
test -f "$prefix/include/seams/seams_lang.h"
test -f "$prefix/lib/cmake/seams/seams-config.cmake"
test -f "$prefix/lib/cmake/seams/seams-config-version.cmake"

for name in seams.pc seams-core.pc seams-lang.pc readcon-core.pc readcon-db.pc; do
  pc=$prefix/lib/pkgconfig/$name
  test -f "$pc"
  line=$(grep '^prefix=' "$pc")
  if [[ $line != 'prefix=${pcfiledir}/../..' ]]; then
    printf '%s\n' "$name prefix is $line" >&2
    exit 1
  fi
done

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
export LD_LIBRARY_PATH=$prefix/lib
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

case $CHECK in
  c)
    if [[ -z ${QUICKSTART:-} || ! -f $QUICKSTART ]]; then
      printf '%s\n' "QUICKSTART must name quickstart.c" >&2
      exit 1
    fi
    cdir=$(mktemp -d)
    # sqrt comes from libm. The seams include and link flags come from pkg-config.
    "$compiler" "$QUICKSTART" -o "$cdir/quickstart" $(pkg-config --cflags --libs seams) -lm
    "$cdir/quickstart"
    printf '%s\n' "c quickstart ok"
    libdir=$(pkg-config --variable=libdir seams)
    "$compiler" "$QUICKSTART" -o "$cdir/quickstart-rpath" \
      $(pkg-config --cflags --libs seams) -lm -Wl,-rpath,"$libdir"
    dyn=$(readelf -d "$cdir/quickstart-rpath")
    printf '%s\n' "$dyn"
    case $dyn in
      *"$libdir"*) ;;
      *)
        printf '%s\n' "quickstart rpath does not name $libdir" >&2
        exit 1
        ;;
    esac
    env -u LD_LIBRARY_PATH "$cdir/quickstart-rpath"
    printf '%s\n' "c quickstart rpath ok"
    ;;
  rust)
    rdir=$(mktemp -d)
    (
      cd "$rdir"
      cargo new readme_example
      cd readme_example
      cargo add dseams@2.14.0
      if grep -q bundled Cargo.toml; then
        printf '%s\n' "dseams was added with the bundled feature" >&2
        exit 1
      fi
    )
    cat > "$rdir/readme_example/src/main.rs" <<'EOF'
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
      cd "$rdir/readme_example"
      cargo build
    )
    example=$rdir/readme_example/target/debug/readme_example
    ldd_out=$(ldd "$example")
    printf '%s\n' "$ldd_out"
    case $ldd_out in
      *"$prefix/lib/libseams.so"*) ;;
      *)
        printf '%s\n' "the example is not linked to the archive libseams" >&2
        exit 1
        ;;
    esac
    "$example"
    printf '%s\n' "dseams readme example ok"
    ;;
  fortran)
    curl -fsSL -o "$stage/fpm.F90" \
      https://github.com/fortran-lang/fpm/releases/download/current/fpm.F90
    gfortran -o "$stage/fpm" "$stage/fpm.F90"
    git clone --depth 1 --branch v2.14.0 \
      https://github.com/d-SEAMS/dseams-fortran.git \
      "$stage/dseams-fortran"
    export FPM_LDFLAGS="$(pkg-config --libs-only-L seams) -Wl,-rpath,$(pkg-config --variable=libdir seams)"
    export OMP_NUM_THREADS=1
    (
      cd "$stage/dseams-fortran"
      "$stage/fpm" test
    )
    printf '%s\n' "dseams-fortran fpm test ok"
    ;;
esac
