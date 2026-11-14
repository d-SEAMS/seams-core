#!/usr/bin/env python3
# Copyright (c) d-SEAMS developers
"""Hold libseams_lang to the C ABI of the 2.13.0 release.

Usage: abi_check_lang.py check|selftest CORPUS SUPPRESSIONS HEADERS LIB WORK

2.13.0 has no libseams; its C interface is bindings/cxx/seams_lang.h in
libseams_lang. CORPUS is that library as libabigail saw it. LIB is the
current libseams_lang built with -g, HEADERS the current seams_lang.h
directory, and SUPPRESSIONS keeps the comparison to the seams_* functions.

check fails when a 2.13.0 function is removed or its type changes. Added
functions pass. selftest compares the corpus with copies of itself, one
with a 2.13.0 function removed and one where a function takes another's
signature, and requires the same abidiff call to fail on each.

SEAMS_ABIDIFF names abidiff. Without it the test is skipped (exit 77), and
so is an abidiff that cannot read the corpus. The corpus is format 2.4,
which libabigail 2.9 and 3.0 read. It was written from a 2.13.0 (7f310225)
build with -g by libabigail 2.9 with
  abidw --headers-dir bindings/cxx --suppressions scripts/abi_seams_lang.supp \
    --drop-private-types --drop-undefined-syms --no-corpus-path \
    --no-comp-dir-path --no-show-locs --no-load-undefined-interfaces \
    -o scripts/abi_2.13.0_seams_lang.abi bindings/libseams_lang.so
and the abi-instr path cut to bindings/cxx/seams_lang.cpp.
"""

import os
import re
import shutil
import subprocess
import sys
from pathlib import Path

SKIP = 77


def tool(env_name, name):
    found = os.environ.get(env_name) or shutil.which(name)
    if not found or not os.access(found, os.X_OK):
        print("{} not found; the 2.13.0 ABI check is skipped".format(name))
        raise SystemExit(SKIP)
    return found


def abidiff(abidiff_bin, corpus, supp, headers, new, out, added):
    cmd = [abidiff_bin, "--drop-private-types", "--no-unreferenced-symbols", "--suppressions", str(supp)]
    if not added:
        cmd.append("--no-added-syms")
    if headers is not None:
        cmd += ["--headers-dir2", str(headers)]
    done = subprocess.run(cmd + [str(corpus), str(new)], capture_output=True, text=True)
    report = done.stdout + done.stderr
    out.write_text(report)
    if done.returncode == 1 and "incompatible format version" in report:
        print(report)
        print("{} cannot read the 2.13.0 corpus; the check is skipped".format(abidiff_bin))
        raise SystemExit(SKIP)
    return done.returncode, report


def check(abidiff_bin, corpus, supp, headers, lib, work):
    status, report = abidiff(abidiff_bin, corpus, supp, headers, lib, work / "abi-2.13.txt", False)
    if status != 0:
        print(report)
        print("libseams_lang breaks the 2.13.0 C ABI (abidiff status {})".format(status), file=sys.stderr)
        return 1
    status, report = abidiff(abidiff_bin, corpus, supp, headers, lib, work / "abi-2.13-added.txt", True)
    if status & 11:
        print(report)
        print("abidiff failed (status {})".format(status), file=sys.stderr)
        return 1
    print(report)
    print("libseams_lang keeps the 2.13.0 C ABI; additions only")
    return 0


def element(xml, tag, name):
    match = re.search(r"\n *<" + tag + r" name='" + name + r"'[^>]*?(?:/>|>.*?</" + tag + ">)", xml, re.S)
    if match is None:
        raise SystemExit("no {} for {} in the dump".format(tag, name))
    return match.group(0)


def selftest(abidiff_bin, corpus, supp, work):
    xml = corpus.read_text()
    if check(abidiff_bin, corpus, supp, None, corpus, work) != 0:
        print("the corpus does not match itself", file=sys.stderr)
        return 1
    rings = element(xml, "function-decl", "seams_frame_rings")
    rings_symbol = element(xml, "elf-symbol", "seams_frame_rings")
    nop = element(xml, "function-decl", "seams_frame_nop")
    retyped = rings.replace("'seams_frame_rings'", "'seams_frame_nop'")
    cases = [
        ("seams_frame_rings removed", xml.replace(rings, "").replace(rings_symbol, "")),
        ("seams_frame_nop with the type of seams_frame_rings", xml.replace(nop, retyped)),
    ]
    failures = []
    for label, text in cases:
        if text == xml:
            failures.append("{}: the edit changed nothing".format(label))
            continue
        tampered = work / "abi-2.13-tampered.abi"
        tampered.write_text(text)
        status, report = abidiff(abidiff_bin, corpus, supp, None, tampered, work / "abi-2.13-tampered.txt", False)
        if status & 12 == 0:
            failures.append("{}: abidiff status {}\n{}".format(label, status, report))
    if failures:
        print("\n".join(failures), file=sys.stderr)
        return 1
    print("the 2.13.0 ABI check rejects a removed function and a changed signature")
    return 0


def main(argv):
    if len(argv) != 7 or argv[1] not in ("check", "selftest"):
        print(__doc__, file=sys.stderr)
        return 2
    mode = argv[1]
    corpus, supp, headers, lib, work = (Path(item) for item in argv[2:])
    work.mkdir(parents=True, exist_ok=True)
    abidiff_bin = tool("SEAMS_ABIDIFF", "abidiff")
    if mode == "selftest":
        return selftest(abidiff_bin, corpus, supp, work)
    if not lib.is_file():
        print("no libseams_lang with debug info at {}; the check is skipped".format(lib))
        return SKIP
    return check(abidiff_bin, corpus, supp, headers, lib, work)


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
