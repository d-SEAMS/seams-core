#!/usr/bin/env python3
# Copyright (c) d-SEAMS developers
"""capi_tsan_gate.py must pass a race inside the main thread's libgomp
team and fail a death, a missing "capi suite ok" and any other report.

Each case is a child that prints what a ThreadSanitizer run would print.
"""

import contextlib
import io
import sys
from pathlib import Path

sys.dont_write_bytecode = True
sys.path.insert(0, str(Path(__file__).resolve().parent))
import capi_tsan_gate  # noqa: E402

BAR = "==================\n"
TEAM_RACE = BAR + """WARNING: ThreadSanitizer: data race (pid=1099568)
  Read of size 8 at 0x7fffffff7898 by thread T1:
    #0 chill::cutoffSolid(...) [clone ._omp_fn.0] ../src/steinhardt_parallel.cpp:1146 (libyodaLib.so+0x4ca2c8)
    #1 gomp_thread_start ../../../libgomp/team.c:129 (libgomp.so.1+0x23f8d)

  Previous write of size 8 at 0x7fffffff7898 by main thread:
    #0 chill::cutoffSolid(...) ../src/steinhardt_parallel.cpp:1143 (libyodaLib.so+0x4e1225)
    #1 main ../tests/capi_suite.c:985 (capi_suite+0x4751)

  Location is stack of main thread.

  Thread T1 (tid=1099570, running) created by main thread at:
    #0 pthread_create ../../../../libsanitizer/tsan/tsan_interceptors_posix.cpp:1022 (libtsan.so.2+0x5baa6)
    #1 gomp_team_start ../../../libgomp/team.c:859 (libgomp.so.1+0x2459b)
    #2 seams_q4q6 ../src/seams_api.cpp:1274 (libseams.so+0x42142)

SUMMARY: ThreadSanitizer: data race ../src/steinhardt_parallel.cpp:1146 in chill::cutoffSolid(...)
""" + BAR
# The same race, where the team belongs to a region another thread hosts.
HOSTED_RACE = TEAM_RACE.replace("created by main thread", "created by thread T4")
THREAD_RACE = BAR + """WARNING: ThreadSanitizer: data race (pid=1177572)
  Write of size 8 at 0x5555555580a0 by main thread:
    #0 main /tmp/gatecheck/stdthread.cpp:6 (stdthread+0x1358)

  Previous write of size 8 at 0x5555555580a0 by thread T1:
    #0 operator() /tmp/gatecheck/stdthread.cpp:5 (stdthread+0x1243)
    #1 execute_native_thread_routine ../../../../../libstdc++-v3/src/c++11/thread.cc:105 (libstdc++.so.6+0xf5b85)

  Location is global 'shared' of size 8 at 0x5555555580a0 (stdthread+0x40a0)

  Thread T1 (tid=1177574, finished) created by main thread at:
    #0 pthread_create ../../../../libsanitizer/tsan/tsan_interceptors_posix.cpp:1022 (libtsan.so.2+0x5baa6)
    #1 std::thread::_M_start_thread(...) gthr-default.h:709 (libstdc++.so.6+0xf5c6e)

SUMMARY: ThreadSanitizer: data race /tmp/gatecheck/stdthread.cpp:6 in main
""" + BAR
USE_AFTER_FREE = BAR + """WARNING: ThreadSanitizer: heap-use-after-free (pid=4242)
  Read of size 8 at 0x7b0400000000 by main thread:
    #0 main ../tests/capi_suite.c:985 (capi_suite+0x4751)

SUMMARY: ThreadSanitizer: heap-use-after-free ../tests/capi_suite.c:985 in main
""" + BAR
DEADLY = """ThreadSanitizer:DEADLYSIGNAL
==4243==ERROR: ThreadSanitizer: SEGV on unknown address 0x000000000000 (pc 0x1 bp 0x2 sp 0x3 T4243)
SUMMARY: ThreadSanitizer: SEGV (capi_cpp_ref+0x115d) in main
==4243==ABORTING
"""
MAPPING = "FATAL: ThreadSanitizer: unexpected memory mapping 0x5f11454b6000-0x5f11454b9000\n"
OK = "capi suite ok"


def child(stdout, stderr, end="raise SystemExit(0)"):
    code = "import sys; sys.stdout.write({!r}); sys.stderr.write({!r}); sys.stdout.flush(); {}".format(
        stdout, stderr, end
    )
    return [sys.executable, "-c", code]


CASES = [
    ("a race inside the main thread's libgomp team", child(OK + "\n", TEAM_RACE * 2), 0),
    ("a clean run", child(OK + "\n", ""), 0),
    ("no sentinel", child("capi fixture ice_ic\n", TEAM_RACE), 1),
    ("a nonzero exit", child(OK + "\n", "", "raise SystemExit(1)"), 1),
    (
        "a SIGABRT after the sentinel",
        child(OK + "\n", "", "import os, signal; os.kill(os.getpid(), signal.SIGABRT)"),
        1,
    ),
    ("a deadly signal in a child, exit 0", child(OK + "\n", DEADLY), 1),
    ("an unexpected mapping in a child, exit 0", child(OK + "\n", MAPPING), 1),
    ("a race with a std::thread", child(OK + "\n", TEAM_RACE + THREAD_RACE), 1),
    ("a race in a team another thread hosts", child(OK + "\n", HOSTED_RACE), 1),
    ("a heap-use-after-free", child(OK + "\n", USE_AFTER_FREE), 1),
]


def main():
    failures = []
    for name, cmd, want in CASES:
        with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()):
            got = capi_tsan_gate.main(["capi_tsan_gate.py"] + cmd)
        if got != want:
            failures.append("{}: gate returned {}, expected {}".format(name, got, want))
    for failure in failures:
        print(failure, file=sys.stderr)
    if not failures:
        print("capi tsan gate: {} cases".format(len(CASES)))
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
