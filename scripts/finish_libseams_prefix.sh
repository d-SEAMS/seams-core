#!/bin/bash
# Copyright (c) d-SEAMS developers
# Make an installed libseams prefix safe to relocate, and record the
# licences of the code compiled into it.
set -euo pipefail

if [[ $# -ne 1 || -z ${EIGEN_SRC:-} ]]; then
  printf '%s\n' "usage: EIGEN_SRC=... finish_libseams_prefix.sh PREFIX" >&2
  exit 1
fi

prefix=$1
root=$(CDPATH= cd -- "$(dirname "$0")/.." && pwd)
shopt -s nullglob
pc_files=("$prefix"/lib/pkgconfig/*.pc)
have_seams=0
for pc in "${pc_files[@]}"; do
  base=$(basename "$pc")
  case $base in
    seams.pc|seams-core.pc|seams-lang.pc) have_seams=$((have_seams + 1)) ;;
    readcon-core.pc|readcon-db.pc) ;;
    *)
      printf '%s\n' "unexpected pkg-config file $base" >&2
      exit 1
      ;;
  esac
  sed -i 's|^prefix=.*|prefix=${pcfiledir}/../..|' "$pc"
  line=$(grep '^prefix=' "$pc")
  if [[ $line != 'prefix=${pcfiledir}/../..' ]]; then
    printf '%s\n' "$base still has prefix line: $line" >&2
    exit 1
  fi
  if grep -F "$prefix" "$pc"; then
    printf '%s\n' "$base still names the build prefix" >&2
    exit 1
  fi
done
if [[ $have_seams -ne 3 ]]; then
  printf '%s\n' "expected seams.pc, seams-core.pc and seams-lang.pc" >&2
  exit 1
fi

cp "$root/LICENSE" "$prefix/LICENSE"
lic=$prefix/licenses
mkdir -p "$lic/argum" "$lic/backward-cpp" "$lic/eigen" "$lic/gcc-runtime"
cp "$root/scripts/licenses/argum.txt" "$lic/argum/LICENSE"
cp "$root/scripts/licenses/backward-cpp.txt" "$lic/backward-cpp/LICENSE"

eigen_copied=0
for file in "$EIGEN_SRC"/COPYING*; do
  if [[ -f $file ]]; then
    cp "$file" "$lic/eigen/$(basename "$file")"
    eigen_copied=1
  fi
done
if [[ $eigen_copied -ne 1 || ! -f $lic/eigen/COPYING.MPL2 || ! -f $lic/eigen/COPYING.BSD ]]; then
  printf '%s\n' "Eigen source is missing COPYING.MPL2 or COPYING.BSD" >&2
  exit 1
fi

copy_gcc_licence() {
  local file=$1
  local base
  base=$(basename "$file")
  case $base in
    COPYING|COPYING3|COPYING.RUNTIME|COPYING3.LIB|COPYING.LIB)
      if [[ -f $file ]]; then
        cp "$file" "$lic/gcc-runtime/$base"
      fi
      ;;
  esac
}

gcc_pkg=$(rpm -qf "$(command -v gcc)")
while IFS= read -r file; do
  copy_gcc_licence "$file"
done < <(rpm -ql "$gcc_pkg")

# manylinux installs that package with documentation omitted, so the
# licence texts are in the rpm payload and absent from the disk.
if [[ ! -f $lic/gcc-runtime/COPYING3 || ! -f $lic/gcc-runtime/COPYING.RUNTIME ]]; then
  if [[ -z $(command -v cpio) || -z $(command -v rpm2cpio) ]]; then
    dnf install -y cpio rpm
  fi
  stage=$(mktemp -d)
  (
    cd "$stage"
    dnf download -y "$gcc_pkg"
    shopt -s nullglob
    rpms=(./*.rpm)
    if [[ ${#rpms[@]} -ne 1 ]]; then
      printf '%s\n' "expected one rpm for $gcc_pkg" >&2
      exit 1
    fi
    rpm2cpio "${rpms[0]}" | cpio -id --quiet
  )
  while IFS= read -r file; do
    copy_gcc_licence "$file"
  done < <(find "$stage" -type f -name 'COPYING*')
  rm -rf "$stage"
fi
if [[ ! -f $lic/gcc-runtime/COPYING3 || ! -f $lic/gcc-runtime/COPYING.RUNTIME ]]; then
  printf '%s\n' "GCC package $gcc_pkg has no runtime licence" >&2
  exit 1
fi

if [[ -n ${READCON_CORE_SRC:-} ]]; then
  mkdir -p "$lic/readcon-core"
  cp "$READCON_CORE_SRC/LICENSE" "$lic/readcon-core/LICENSE"
fi
if [[ -n ${READCON_DB_SRC:-} ]]; then
  mkdir -p "$lic/readcon-db"
  cp "$READCON_DB_SRC/LICENSE" "$lic/readcon-db/LICENSE"
fi
if [[ ! -f $lic/readcon-core/LICENSE || ! -f $lic/readcon-db/LICENSE ]]; then
  printf '%s\n' "readcon licence texts were not copied" >&2
  exit 1
fi

cat > "$lic/README" <<'EOF'
The project licence is LICENSE at the top of this archive.

These directories hold the licences of code compiled into the binaries.

eigen
  Eigen headers compiled into libyodaLib. MPL-2.0, and the BSD terms
  in COPYING.BSD for the files that document names.

argum
  The argument parser compiled into the seams command. BSD-3-Clause.

backward-cpp
  The crash handler compiled into the seams command. MIT.

gcc-runtime
  GCC runtime objects the manylinux_2_28 compiler links statically
  for symbols absent from the shared libstdc++ and libgcc on glibc
  2.28. GPL-3.0 with the GCC Runtime Library Exception.

readcon-core
  The .con reader shipped beside libseams. MIT.

readcon-db
  The CON corpus reader shipped beside libseams. MIT.
EOF
