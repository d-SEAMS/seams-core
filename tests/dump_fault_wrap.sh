#!/bin/sh
# Copyright (c) d-SEAMS developers
# The fault library has to be preloaded, and libasan has to stay ahead of
# it. Meson's test environment would replace LD_PRELOAD with only the
# fault library, and the driver then aborts before main.
fault=$1
shift
if [ -n "${ASAN_TEST_PRELOAD:-}" ]; then
  export LD_PRELOAD="${ASAN_TEST_PRELOAD}:$fault"
else
  export LD_PRELOAD="${LD_PRELOAD:+$LD_PRELOAD:}$fault"
fi
exec "$@"
