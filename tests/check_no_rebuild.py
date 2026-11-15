#!/usr/bin/env python3
# Copyright (c) d-SEAMS developers
"""No test rebuilds the tree that the other tests are loading.

meson test runs tests in parallel against one build directory. A test
that calls ``meson install`` without ``--no-rebuild``, ``meson compile``
or ``ninja`` on that directory can relink libyodaLib.so or libseams.so
while another test has it mapped. This check fails when a test script
does that, and when c_api_installed installs a target that meson test
would not have built first. A test that reads BUILD/docs-r-lib has to
depend on the docs-r-lib target, so no test installs the R package
while another test loads it.

usage: check_no_rebuild.py SOURCE_ROOT BUILD_DIR
"""

import json
import re
import shutil
import subprocess
import sys
from pathlib import Path

CALL = re.compile(r"""["']?\bmeson["']?\s*,?\s*["']?(install|compile)\b|["']ninja["']|\bninja\s+-C\b""")
INSTALL_TEST = "c_api_installed"
R_LIB = "docs-r-lib"


def call_problems(text):
    """Rebuilding calls in a script body: (line number, line) pairs."""
    found = []
    lines = text.splitlines()
    for number, line in enumerate(lines, start=1):
        stripped = line.lstrip()
        if stripped.startswith("#"):
            continue
        for match in CALL.finditer(line):
            verb = match.group(1)
            if verb == "install":
                call = " ".join(lines[number - 1:number + 2])
                if "--no-rebuild" in call:
                    continue
            found.append((number, line.strip()))
            break
    return found


def selftest():
    bad = [
        'run(["meson", "install", "-C", str(build), "--destdir", str(stage)])',
        'subprocess.run(["meson", "compile", "-C", build])',
        'subprocess.run(["ninja", "-C", build, "src/libyodaLib.so"])',
        "ninja -C build src/libseams.so",
    ]
    good = [
        'run(["meson", "install", "-C", str(build), "--no-rebuild", "--destdir", str(stage)])',
        "# meson compile -C build is how the tree was built",
        '["meson", "introspect", "--buildoptions", str(build)]',
    ]
    for text in bad:
        if not call_problems(text):
            raise SystemExit(f"selftest: missed {text}")
    for text in good:
        if call_problems(text):
            raise SystemExit(f"selftest: flagged {text}")


def introspect(build):
    meson = shutil.which("meson") or "meson"
    out = subprocess.check_output([meson, "introspect", "--targets", "--tests", str(build)], text=True)
    data = json.loads(out)
    return data["targets"], data["tests"]


def test_scripts(tests, source):
    scripts = {}
    for test in tests:
        for item in list(test["cmd"]) + list(test.get("args", [])):
            path = Path(item)
            if path.suffix not in (".py", ".sh"):
                continue
            try:
                path.resolve().relative_to(source)
            except ValueError:
                continue
            if path.is_file():
                scripts.setdefault(path.resolve(), []).append(test["name"])
    return scripts


def ninja_inputs(build, outputs):
    ninja = shutil.which("ninja")
    if ninja is None or not outputs:
        return set()
    done = subprocess.run(
        [ninja, "-C", str(build), "-t", "inputs", *outputs],
        text=True,
        capture_output=True,
    )
    if done.returncode != 0:
        return set()
    return {str((build / line.strip()).resolve()) for line in done.stdout.splitlines() if line.strip()}


def install_problems(build, targets, tests):
    rows = [test for test in tests if test["name"] == INSTALL_TEST]
    if not rows:
        return []
    by_id = {target["id"]: target for target in targets}
    built = set()
    outputs = []
    for target_id in rows[0]["depends"]:
        target = by_id.get(target_id)
        if target is None:
            continue
        for name in target["filename"]:
            built.add(str(Path(name).resolve()))
            outputs.append(str(Path(name).resolve().relative_to(build.resolve())))
    built |= ninja_inputs(build.resolve(), outputs)
    missing = []
    for target in targets:
        if not target.get("installed"):
            continue
        for name in target["filename"]:
            if str(Path(name).resolve()) not in built:
                missing.append(target["name"])
                break
    return sorted(set(missing))


def r_lib_problems(targets, tests, source):
    """Tests whose scripts read docs-r-lib without depending on its target."""
    ids = [target["id"] for target in targets if target["name"] == R_LIB]
    if not ids:
        return []
    missing = []
    for test in tests:
        if ids[0] in test["depends"]:
            continue
        for script in test_scripts([test], source):
            if script == Path(__file__).resolve():
                continue
            if R_LIB in script.read_text(encoding="utf-8", errors="replace"):
                missing.append(test["name"])
                break
    return sorted(set(missing))


def main():
    if len(sys.argv) != 3:
        raise SystemExit(__doc__)
    selftest()
    source = Path(sys.argv[1]).resolve()
    build = Path(sys.argv[2])
    targets, tests = introspect(build)
    failures = []
    for script, names in sorted(test_scripts(tests, source).items()):
        if script == Path(__file__).resolve():
            continue
        for number, line in call_problems(script.read_text(encoding="utf-8", errors="replace")):
            rel = script.relative_to(source)
            failures.append(f"{rel}:{number} rebuilds the tree under {', '.join(sorted(set(names)))}: {line}")
    for name in install_problems(build, targets, tests):
        failures.append(f"{INSTALL_TEST} installs {name} but does not depend on it")
    for name in r_lib_problems(targets, tests, source):
        failures.append(f"{name} reads {R_LIB} but does not depend on it")
    if failures:
        print("\n".join(failures))
        return 1
    print(f"no test rebuilds the tree; {INSTALL_TEST} depends on every installed target")
    return 0


if __name__ == "__main__":
    sys.exit(main())
