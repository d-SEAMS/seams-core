#!/usr/bin/env python3
# Copyright (c) d-SEAMS developers
"""ylmAdd8 must not be referenced without the AVX-512 object.

chill::cutoffSolid calls seams::ylm::ylmAdd8. That function is defined
only in ylm_cartesian_avx512.cpp, and meson builds that file only on
x86-64. The call used to sit behind a flag that is constant false on
every other machine. Apple clang at -O0 keeps the call, so a debug
arm64 link fails with an undefined symbol. The archive named here is
steinhardt_parallel.cpp compiled at -O0 with SEAMS_YLM_AVX512 unset.
"""

import subprocess
import sys


def main():
    if len(sys.argv) != 2:
        sys.exit("usage: check_ylm_add8.py ARCHIVE")
    done = subprocess.run(
        ["nm", "-u", sys.argv[1]],
        check=False,
        capture_output=True,
        text=True,
    )
    if done.returncode != 0:
        sys.stderr.write(done.stderr)
        sys.exit("nm failed on {}".format(sys.argv[1]))
    hits = [line for line in done.stdout.splitlines() if "ylmAdd8" in line]
    if hits:
        sys.stderr.write("\n".join(hits) + "\n")
        sys.exit(
            "ylmAdd8 is undefined at -O0 without SEAMS_YLM_AVX512"
        )
    print("ylmAdd8 is not referenced at -O0 without AVX-512")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
