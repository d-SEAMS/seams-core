#!/bin/sh
# Copyright (c) d-SEAMS developers
# Run one cargo test binary with the sanitizer runtime already mapped.
# rustc must not see this preload: cc1 then leak-checks itself and the
# cargo-test budget expires during the compile. The parent cargo wrapper
# drops LD_PRELOAD and points CARGO_TARGET_*_RUNNER here, and only when
# the test loads libseams.
if [ -n "${ASAN_TEST_PRELOAD:-}" ]; then
  case ":${LD_PRELOAD:-}:" in
    *":${ASAN_TEST_PRELOAD}:"*) ;;
    *) LD_PRELOAD="${ASAN_TEST_PRELOAD}${LD_PRELOAD:+:$LD_PRELOAD}" ;;
  esac
  export LD_PRELOAD
fi
exec "$@"
