#!/usr/bin/env -S -u LD_PRELOAD /bin/sh
# Run one cargo test binary with the sanitizer runtime already mapped.
# rustc must not see this preload: cc1 then leak-checks itself and the
# cargo-test budget expires during the compile. The parent cargo wrapper
# drops LD_PRELOAD and points CARGO_TARGET_*_RUNNER here, and only when
# the test loads libseams.
# cargo itself must not be preloaded (rustc and /bin/sh die). This
# runner is exec'd only for the test or example. A Rust binary that
# links libseams does not list libclang_rt.tsan in its own NEEDED, so
# the runtime has to be here, before that process starts. A binary that
# already needs the shared runtime still accepts this preload: it is
# the same file.
if [ -n "${TSAN_RUNTIME:-}" ] && [ -f "$TSAN_RUNTIME" ]; then
  case ":${LD_PRELOAD:-}:" in
    *":${TSAN_RUNTIME}:"*) ;;
    *) LD_PRELOAD="${TSAN_RUNTIME}${LD_PRELOAD:+:$LD_PRELOAD}" ;;
  esac
  export LD_PRELOAD
fi
if [ -n "${ASAN_TEST_PRELOAD:-}" ]; then
  case ":${LD_PRELOAD:-}:" in
    *":${ASAN_TEST_PRELOAD}:"*) ;;
    *) LD_PRELOAD="${ASAN_TEST_PRELOAD}${LD_PRELOAD:+:$LD_PRELOAD}" ;;
  esac
  export LD_PRELOAD
fi
exec "$@"
