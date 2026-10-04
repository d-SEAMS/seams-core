#!/bin/sh
# Profile sampleRDF_AA with the POP efficiencies.
# Load balance, communication efficiency, and parallel efficiency come from
# a GOMP_parallel interposer. IPC is reported when the PMU opens.
set -eu
root=$(CDPATH= cd -- "$(dirname "$0")/../.." && pwd)
cc=${CC:-/workspace/.pixi/envs/default/bin/x86_64-conda-linux-gnu-gcc}
cxx=${CXX:-/workspace/.pixi/envs/default/bin/x86_64-conda-linux-gnu-c++}
if [ ! -x "$cc" ]; then cc=gcc; fi
if [ ! -x "$cxx" ]; then cxx=c++; fi
pixi_prefix=$(CDPATH= cd -- "$(dirname "$cc")/.." && pwd)
out=${POP3_OUT:-/tmp/pop3-rdf}
mkdir -p "$out"
"$cc" -shared -fPIC -O2 -fopenmp \
  -I"$root/analysis/pop3" \
  "$root/analysis/pop3/gomp_pop.c" \
  -Wl,--version-script="$root/analysis/pop3/gomp_pop.map" \
  -ldl -lgomp \
  -o "$out/libgomp_pop.so"
"$cxx" -std=c++20 -O2 -fopenmp \
  -I"$root/analysis/pop3" \
  -I"$root/src/include/internal" \
  -I"$root/src/include/external" \
  -I"$root/subprojects/minimage/include" \
  -I"$pixi_prefix/include/eigen3" \
  -DSEAMS_HAS_OPENMP=1 \
  -DSEAMS_HAS_MINIMAGE=1 \
  "$root/analysis/pop3/rdf_pop.cpp" \
  -L"$root/bbdir/src" -lyodaLib \
  -Wl,-rpath,"$root/bbdir/src:$root/bbdir/subprojects/readcon-db:$pixi_prefix/lib" \
  -Wl,-rpath-link,"$root/bbdir/subprojects/readcon-db:$root/bbdir/subprojects/readcon-core" \
  -ldl \
  -o "$out/rdf_pop"
LD_PRELOAD="$out/libgomp_pop.so" "$out/rdf_pop"
