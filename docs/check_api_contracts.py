#!/usr/bin/env python3
# Copyright (c) d-SEAMS developers
"""Verify capi_docs.py --write and write_api_contracts.py are fixed points.

Runs both generators on temporary copies of docs/api/c.md and
docs/api/cpp.md and fails on any byte difference against the committed
pages, then exercises negative self-test mutations to confirm that drift
on either page is detected and repaired.
"""

import difflib
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DOCS = ROOT / "docs"
API = DOCS / "api"
HEADER = ROOT / "include" / "seams.h"
SYMBOLS = ROOT / "tests" / "capi_symbols.txt"
CAPI_DOCS = DOCS / "capi_docs.py"
WRITE_CONTRACTS = DOCS / "write_api_contracts.py"
PAGES = ("c.md", "cpp.md")


def _run(cmd):
    return subprocess.run(cmd, cwd=str(ROOT), check=False, capture_output=True, text=True)


def _run_generators(work_dir):
    step1 = _run(
        [
            sys.executable,
            str(CAPI_DOCS),
            "--write",
            str(HEADER),
            str(SYMBOLS),
            str(work_dir / "c.md"),
        ]
    )
    if step1.returncode != 0:
        return "capi_docs.py --write failed ({}):\n{}{}".format(
            step1.returncode, step1.stdout, step1.stderr
        )
    step2 = _run(
        [
            sys.executable,
            str(WRITE_CONTRACTS),
            "--api-dir",
            str(work_dir),
        ]
    )
    if step2.returncode != 0:
        return "write_api_contracts.py failed ({}):\n{}{}".format(
            step2.returncode, step2.stdout, step2.stderr
        )
    return None


def _copy_pages(work_dir):
    for name in PAGES:
        shutil.copy2(API / name, work_dir / name)


def _unified_diff(name, before, after):
    return "".join(
        difflib.unified_diff(
            before.splitlines(keepends=True),
            after.splitlines(keepends=True),
            fromfile="docs/api/" + name,
            tofile="generated/" + name,
        )
    )


def _check_fixed_point(failures):
    with tempfile.TemporaryDirectory(prefix="seams-api-contracts-") as tmp:
        work = Path(tmp)
        _copy_pages(work)
        err = _run_generators(work)
        if err:
            failures.append(err)
            return
        for name in PAGES:
            committed = (API / name).read_text()
            generated = (work / name).read_text()
            if generated != committed:
                failures.append(
                    "docs/api/{} changed after running generators:\n{}".format(
                        name, _unified_diff(name, committed, generated)
                    )
                )

    check = _run([sys.executable, str(WRITE_CONTRACTS), "--check"])
    if check.returncode != 0:
        failures.append(
            "write_api_contracts.py --check failed ({}):\n{}{}".format(
                check.returncode, check.stdout, check.stderr
            )
        )


def _selftest(failures):
    committed_c = (API / "c.md").read_text()
    committed_cpp = (API / "cpp.md").read_text()

    # 1. Missing Errors line in c.md must be detected and restored.
    with tempfile.TemporaryDirectory(prefix="seams-api-selftest-c-") as tmp:
        work = Path(tmp)
        _copy_pages(work)
        mutated = re.sub(
            r"(### seams_tetrahedral_order\n.*?)(\n\nErrors: [^\n]+)",
            r"\1",
            committed_c,
            count=1,
            flags=re.S,
        )
        if mutated == committed_c:
            failures.append("selftest(c): failed to mutate seams_tetrahedral_order")
        else:
            (work / "c.md").write_text(mutated)
            chk = _run([sys.executable, str(WRITE_CONTRACTS), "--check", "--api-dir", str(work)])
            if chk.returncode == 0:
                failures.append("selftest(c): --check missed a dropped Errors line in c.md")
            err = _run_generators(work)
            if err:
                failures.append("selftest(c): " + err)
            elif (work / "c.md").read_text() != committed_c:
                failures.append("selftest(c): generators did not restore committed c.md")

    # 2. Stale Purpose text in c.md must be rewritten back from seams.h by capi_docs.py --write.
    with tempfile.TemporaryDirectory(prefix="seams-api-selftest-c-hdr-") as tmp:
        work = Path(tmp)
        _copy_pages(work)
        mutated = re.sub(
            r"(### seams_tetrahedral_order\n\nPurpose: )([^\n]+)",
            r"\1Stale purpose text.",
            committed_c,
            count=1,
        )
        if mutated == committed_c:
            failures.append("selftest(c-hdr): failed to mutate seams_tetrahedral_order purpose")
        else:
            (work / "c.md").write_text(mutated)
            err = _run_generators(work)
            if err:
                failures.append("selftest(c-hdr): " + err)
            else:
                repaired = (work / "c.md").read_text()
                if repaired == mutated or repaired != committed_c:
                    failures.append("selftest(c-hdr): capi_docs.py --write did not restore committed c.md")

    # 3. Missing Returns line on seams_frame_read_options in cpp.md must be detected and filled.
    with tempfile.TemporaryDirectory(prefix="seams-api-selftest-cpp-ret-") as tmp:
        work = Path(tmp)
        _copy_pages(work)
        mutated = re.sub(
            r"(### seams_frame_read_options\n.*?)(\n\nReturns: [^\n]+)",
            r"\1",
            committed_cpp,
            count=1,
            flags=re.S,
        )
        if mutated == committed_cpp:
            failures.append("selftest(cpp-ret): failed to mutate seams_frame_read_options")
        else:
            (work / "cpp.md").write_text(mutated)
            chk = _run([sys.executable, str(WRITE_CONTRACTS), "--check", "--api-dir", str(work)])
            if chk.returncode == 0:
                failures.append("selftest(cpp-ret): --check missed a dropped Returns line in cpp.md")
            err = _run_generators(work)
            if err:
                failures.append("selftest(cpp-ret): " + err)
            else:
                repaired = (work / "cpp.md").read_text()
                sec = re.search(
                    r"### seams_frame_read_options\n.*?(?=\n### |\n<!-- API_SYMBOLS_END -->)",
                    repaired,
                    re.S,
                )
                expected = "Returns: A new seams_frame, or NULL when the call cannot build one."
                if not sec or expected not in sec.group(0):
                    failures.append(
                        "selftest(cpp-ret): missing seams_frame * return contract for seams_frame_read_options"
                    )

    # 4. Missing symbol section seams_graph_rings in cpp.md must be detected and regenerated.
    with tempfile.TemporaryDirectory(prefix="seams-api-selftest-cpp-miss-") as tmp:
        work = Path(tmp)
        _copy_pages(work)
        mutated = re.sub(
            r"\n### seams_graph_rings\n.*?(?=\n<!-- API_SYMBOLS_END -->)",
            "",
            committed_cpp,
            count=1,
            flags=re.S,
        )
        if mutated == committed_cpp:
            failures.append("selftest(cpp-miss): failed to remove seams_graph_rings")
        else:
            (work / "cpp.md").write_text(mutated)
            chk = _run([sys.executable, str(WRITE_CONTRACTS), "--check", "--api-dir", str(work)])
            if chk.returncode == 0:
                failures.append("selftest(cpp-miss): --check missed removed seams_graph_rings in cpp.md")
            err = _run_generators(work)
            if err:
                failures.append("selftest(cpp-miss): " + err)
            else:
                repaired = (work / "cpp.md").read_text()
                sec = re.search(
                    r"### seams_graph_rings\n.*?(?=\n### |\n<!-- API_SYMBOLS_END -->)",
                    repaired,
                    re.S,
                )
                expected_own = "Ownership: The caller owns every input and output buffer."
                if not sec or expected_own not in sec.group(0):
                    failures.append(
                        "selftest(cpp-miss): regenerated seams_graph_rings lacks non-frame ownership"
                    )

    # 5. Obsolete symbol section in cpp.md must be detected and pruned to exact committed cpp.md.
    with tempfile.TemporaryDirectory(prefix="seams-api-selftest-cpp-stale-") as tmp:
        work = Path(tmp)
        _copy_pages(work)
        stale = (
            "\n### seams_nonexistent_symbol\n\n"
            "Purpose: Obsolete entry.\n\n"
            "Parameters: none\n\n"
            "Returns: None.\n\n"
            "Ownership: No allocation.\n\n"
            "Errors: None.\n\n"
            "Units: dimensionless\n\n"
            "Example: `docs/examples/minimal_cpp.cpp`\n"
        )
        mutated = committed_cpp.replace("<!-- API_SYMBOLS_END -->", stale + "<!-- API_SYMBOLS_END -->", 1)
        (work / "cpp.md").write_text(mutated)
        chk = _run([sys.executable, str(WRITE_CONTRACTS), "--check", "--api-dir", str(work)])
        if chk.returncode == 0:
            failures.append("selftest(cpp-stale): --check missed obsolete symbol in cpp.md")
        err = _run_generators(work)
        if err:
            failures.append("selftest(cpp-stale): " + err)
        elif (work / "cpp.md").read_text() != committed_cpp:
            failures.append("selftest(cpp-stale): pruning obsolete symbol did not reproduce committed cpp.md")


def main():
    failures = []
    _check_fixed_point(failures)
    _selftest(failures)
    if failures:
        print("\n".join(failures), file=sys.stderr)
        return 1
    print("api contracts fixed point ok (c.md, cpp.md, 5 selftest cases)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
