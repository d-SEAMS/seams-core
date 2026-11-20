#!/usr/bin/env python3
# Copyright (c) d-SEAMS developers
"""Check that configure never keeps a chemfiles with fewer molfile readers.

    test_chemfiles_configure.py MESON_LOG LINKED_PROBE [SYSTEM_PROBE]

LINKED_PROBE is cmake/chemfiles_mmcif_probe.cpp built against the chemfiles
this build links; it must pass. SYSTEM_PROBE is the same source built
against the installed chemfiles configure found, when it found one. When
that library lacks molfile readers (probe exit 7), configure must have
printed how many it has and built the patched copy. When it passes,
configure must have kept it.
"""
import re
import subprocess
import sys

EXPECTED = 33


def run(probe):
    return subprocess.run([probe], capture_output=True, text=True, timeout=120)


def main():
    if len(sys.argv) not in (3, 4):
        sys.exit(__doc__)
    log = open(sys.argv[1], encoding="utf-8", errors="replace").read()
    # A reconfigure appends to the log; only the last configure counts.
    log = log.split("Build started at")[-1]

    linked = run(sys.argv[2])
    if linked.returncode != 0:
        sys.exit(
            "the linked chemfiles fails the probe (exit {}): {}{}".format(
                linked.returncode, linked.stdout, linked.stderr))

    if len(sys.argv) == 3:
        if "Message: system chemfiles" in log:
            sys.exit("configure judged a system chemfiles but none was passed here")
        print("no system chemfiles; the linked one has all {} molfile readers".format(EXPECTED))
        return

    system = run(sys.argv[3])
    if system.returncode == 7:
        found = re.fullmatch(
            r"(\d+) of {} VMD molfile readers; missing (.+)".format(EXPECTED),
            system.stdout.strip())
        if found is None or int(found.group(1)) >= EXPECTED:
            sys.exit("unexpected probe output: " + system.stdout)
        message = ("Message: system chemfiles has {}; building patched 0.10.4, "
                   "which has all of them").format(system.stdout.strip())
        if message not in log:
            sys.exit("configure did not report the missing readers:\n" + message)
        print("system chemfiles has {} of {} molfile readers; configure said so "
              "and built the patched copy".format(found.group(1), EXPECTED))
        return
    if system.returncode == 0:
        if "Message: system chemfiles" in log:
            sys.exit("the system chemfiles passes the probe but configure rejected it")
        print("system chemfiles passes the probe and configure kept it")
        return
    if "Message: system chemfiles" not in log or "building patched 0.10.4" not in log:
        sys.exit("system chemfiles fails the probe (exit {}) and configure did not "
                 "say it built the patched copy".format(system.returncode))
    print("system chemfiles fails the probe (exit {}); configure built the patched "
          "copy".format(system.returncode))


if __name__ == "__main__":
    main()
