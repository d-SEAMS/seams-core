#!/usr/bin/env python3
# Copyright (c) d-SEAMS developers
"""Parse every ``.. mermaid::`` block with the mermaid version in conf.py."""

import json
import os
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DOCS = ROOT / "docs"
CONF = DOCS / "source" / "conf.py"
SCRIPT = DOCS / "check_mermaid.mjs"
JSDOM = "25.0.1"
SKIP = {"build", "doctrees", ".cache", "node_modules", "__pycache__"}
DIRECTIVE = re.compile(r"^.. mermaid::\s*$")
OPTION = re.compile(r"^:[A-Za-z0-9_-]+:")


def fail(message):
    print(message, file=sys.stderr)
    raise SystemExit(1)


def mermaid_version():
    match = re.search(
        r'^mermaid_version\s*=\s*"([^"]+)"',
        CONF.read_text(encoding="utf-8"),
        re.M,
    )
    if match is None:
        fail("docs/source/conf.py has no mermaid_version")
    return match.group(1)


def directive_body(lines, start):
    """Return the diagram text and the index of the next line."""
    index = start
    collected = []
    while index < len(lines):
        line = lines[index]
        if line.strip() == "":
            if collected and (index + 1 >= len(lines) or not lines[index + 1].startswith(" ")):
                break
            collected.append("")
            index += 1
            continue
        if not line.startswith(" "):
            break
        collected.append(line)
        index += 1
    if not collected:
        return "", index
    indent = min(len(line) - len(line.lstrip(" ")) for line in collected if line.strip())
    body = []
    for line in collected:
        if line.strip() == "":
            body.append("")
        else:
            body.append(line[indent:])
    while body and body[0] == "":
        body.pop(0)
    while body and OPTION.match(body[0]):
        body.pop(0)
        if body and body[0] == "":
            body.pop(0)
    return "\n".join(body).strip() + "\n", index


def blocks():
    found = []
    for path in DOCS.rglob("*"):
        if not path.is_file() or path.suffix not in {".org", ".rst", ".md"}:
            continue
        if SKIP.intersection(path.relative_to(DOCS).parts):
            continue
        lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
        index = 0
        while index < len(lines):
            if DIRECTIVE.match(lines[index].strip()):
                text, index = directive_body(lines, index + 1)
                if not text.strip():
                    fail("{} has an empty mermaid block".format(path.relative_to(ROOT)))
                found.append({"path": str(path.relative_to(ROOT)), "text": text})
            else:
                index += 1
    return found


def install(version, dest):
    npm = shutil.which("npm")
    if npm is None:
        fail("npm is required to install mermaid {}".format(version))
    done = subprocess.run(
        [
            npm,
            "install",
            "--prefix",
            str(dest),
            "--silent",
            "mermaid@{}".format(version),
            "jsdom@{}".format(JSDOM),
        ],
        cwd=str(ROOT),
        text=True,
        capture_output=True,
    )
    if done.returncode != 0:
        detail = (done.stderr or done.stdout or "").strip()
        fail("npm install mermaid@{} failed\n{}".format(version, detail[-800:]))
    package = dest / "node_modules" / "mermaid"
    jsdom = dest / "node_modules" / "jsdom" / "lib" / "api.js"
    if not (package / "package.json").is_file() or not jsdom.is_file():
        fail("mermaid or jsdom did not install under {}".format(dest))
    return package, jsdom


def main():
    version = mermaid_version()
    found = blocks()
    if not found:
        fail("no .. mermaid:: blocks under docs/")
    with tempfile.TemporaryDirectory(prefix="seams-mermaid-") as temp:
        package, jsdom = install(version, Path(temp))
        env = os.environ.copy()
        env["MERMAID_PACKAGE"] = str(package)
        env["JSDOM_MODULE"] = str(jsdom)
        env["MERMAID_EXPECT"] = version
        done = subprocess.run(
            ["node", str(SCRIPT)],
            input=json.dumps(found),
            text=True,
            capture_output=True,
            env=env,
            cwd=str(ROOT),
        )
    sys.stdout.write(done.stdout)
    sys.stderr.write(done.stderr)
    return done.returncode


if __name__ == "__main__":
    raise SystemExit(main())
