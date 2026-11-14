#!/usr/bin/env python3
# Copyright (c) d-SEAMS developers
"""Execute every Python block and every runnable seams CLI example.

Prints, per page, how many blocks were found, how many ran, and how many
passed. A doctest checks the printed value. A code block with an assert
checks that assert. A seams line whose following text quotes an output
fragment must reproduce that fragment.
"""

from __future__ import annotations

import os
import re
import shutil
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "docs"))
import run_snippets  # noqa: E402

PY_BLOCK = re.compile(
    r"```python\n(.*?)```|#\+begin_src python[^\n]*\n(.*?)#\+end_src|\.\.\s+doctest::\n((?:[ \t][^\n]*\n|\n)*)",
    re.S,
)
BASH_BLOCK = re.compile(
    r"```(?:bash|sh|shell|console)\n(.*?)```|#\+begin_src (?:bash|sh|shell)[^\n]*\n(.*?)#\+end_src",
    re.S,
)
QUOTE = re.compile(r"~([^~\n]+)~|`([^`\n]+)`")
SKIP_LINE = re.compile(r"\b(pip|conda|meson|pixi|nix|mpiexec|EB=|cc |g\+\+|cargo|R CMD)\b")


def pages():
    found = []
    for path in (ROOT / "docs").rglob("*"):
        if path.suffix not in {".org", ".md", ".rst"}:
            continue
        rel = path.relative_to(ROOT / "docs")
        if rel.parts[:1] == ("source",) and (len(rel.parts) < 2 or rel.parts[1] != "lang"):
            continue
        if "videos" in rel.parts:
            continue
        found.append(path)
    return sorted(found)


def python_blocks(text):
    blocks = []
    for match in PY_BLOCK.finditer(text):
        body = next(part for part in match.groups() if part)
        body = run_snippets.unindent(body)
        if body.strip():
            blocks.append(body)
    return blocks


def bash_blocks(text):
    blocks = []
    for match in BASH_BLOCK.finditer(text):
        body = next(part for part in match.groups() if part)
        if re.search(r"(?m)^\s*seams(\s|$)", body):
            blocks.append(body)
    return blocks


def logical_lines(body):
    lines = []
    buf = ""
    for raw in body.splitlines():
        stripped = raw.strip()
        if not stripped or stripped.startswith("#"):
            if buf:
                lines.append(buf)
                buf = ""
            continue
        buf = (buf + " " + stripped).strip() if buf else stripped
        if buf.endswith("\\"):
            buf = buf[:-1].rstrip()
            continue
        lines.append(buf)
        buf = ""
    if buf:
        lines.append(buf)
    return lines


def runnable_seams(body):
    kept = []
    dump = None
    for line in logical_lines(body):
        if line.startswith("dump="):
            dump = line.split("=", 1)[1].strip().strip("\"'")
            continue
        if not line.startswith("seams"):
            continue
        if dump:
            line = line.replace('"$dump"', dump).replace("$dump", dump)
        if SKIP_LINE.search(line) or re.search(r"\bFILE\b|water\.lammpstrj|dump\.lammpstrj", line):
            continue
        missing = False
        for token in line.split():
            token = token.strip("\"'")
            if "/" in token or token.endswith((".lammpstrj", ".xyz", ".gro", ".xtc", ".env")):
                if not (ROOT / token).exists():
                    missing = True
        if missing or line == "seams":
            continue
        kept.append(line)
    return "\n".join(kept) if kept else None


def expected_phrases(text, start):
    window = text[start:start + 900]
    stop = len(window)
    for marker in ("#+begin_src", "```"):
        at = window.find(marker)
        if at > 0:
            stop = min(stop, at)
    phrases = []
    for match in QUOTE.finditer(window[:stop]):
        phrase = re.sub(r"\s+", " ", (match.group(1) or match.group(2) or "").strip())
        if len(phrase) > 80 or "..." in phrase or phrase.startswith("--") or "lx" in phrase:
            continue
        pairs = re.findall(
            r"\b(nop|hexagonal|cubic|water|interClathrate|clathrate|interfacial|rings)\s+(\d+)",
            phrase,
        )
        if not pairs:
            continue
        phrases.append(pairs)
    return phrases


def run_python_block(build, env, body):
    if "```" in body:
        raise run_snippets.Fail("code fence leaked into a Python block")
    if body.lstrip().startswith(">>>") or "\n>>> " in body:
        run_snippets.run_doctest(env, body)
        return
    if body.lstrip().startswith("# excerpt:"):
        target = run_snippets.check_excerpt(body)
        if target is None:
            return
        try:
            if "tests" not in target.parts:
                run_snippets.run_source(build, env, target)
        except run_snippets.Skip:
            return
        return
    run_snippets.run_python(build, env, [body])


def run_cli_block(build, script):
    env = os.environ.copy()
    env["PATH"] = str(build) + os.pathsep + env.get("PATH", "")
    lib = str(build / "src")
    env["LD_LIBRARY_PATH"] = lib + os.pathsep + env.get("LD_LIBRARY_PATH", "")
    done = subprocess.run(
        ["bash", "-c", script],
        cwd=ROOT,
        env=env,
        text=True,
        capture_output=True,
        timeout=180,
    )
    if done.returncode != 0:
        raise run_snippets.Fail((done.stderr or done.stdout)[-500:])
    return done.stdout + done.stderr


def main():
    if len(sys.argv) != 2:
        print("usage: run_doc_blocks.py BUILD", file=sys.stderr)
        return 2
    build = Path(sys.argv[1]).resolve()
    if not (build / "seams").is_file():
        print("seams was not built", file=sys.stderr)
        return 1
    env = run_snippets.base_env(build)
    failures = []
    print(f"{'page':<48} {'found':>6} {'ran':>6} {'pass':>6}")
    for path in pages():
        text = path.read_text(encoding="utf-8")
        py = python_blocks(text)
        cli = bash_blocks(text)
        if not py and not cli:
            continue
        found = len(py) + len(cli)
        ran = 0
        passed = 0
        rel = str(path.relative_to(ROOT))
        for body in py:
            ran += 1
            try:
                run_python_block(build, env, body)
            except run_snippets.Skip as exc:
                failures.append(rel + " python skipped: " + str(exc))
            except Exception as exc:
                failures.append(rel + " python: " + str(exc).splitlines()[-1][:240])
            else:
                passed += 1
        for body in cli:
            script = runnable_seams(body)
            if script is None:
                continue
            ran += 1
            start = text.find(body)
            try:
                output = run_cli_block(build, script)
                for pairs in expected_phrases(text, start + len(body) if start >= 0 else 0):
                    for name, number in pairs:
                        token = name + " " + number
                        if token not in output:
                            raise run_snippets.Fail("missing output " + token)
            except Exception as exc:
                failures.append(rel + " cli: " + str(exc).splitlines()[-1][:240])
            else:
                passed += 1
        print(f"{rel:<48} {found:6d} {ran:6d} {passed:6d}")
    created = ROOT / "seams.env"
    if created.is_file() and "seams.env" not in subprocess.check_output(
        ["git", "-C", str(ROOT), "ls-files", "seams.env"], text=True
    ):
        created.unlink()
    output_dir = ROOT / "output"
    if output_dir.is_dir() and not any(output_dir.iterdir()):
        output_dir.rmdir()
    elif output_dir.is_dir():
        shutil.rmtree(output_dir)
    if failures:
        print("FAILED", len(failures))
        for item in failures:
            print(item)
        return 1
    print("all documented Python blocks and seams examples passed")
    return 0


if __name__ == "__main__":
    sys.exit(main())
