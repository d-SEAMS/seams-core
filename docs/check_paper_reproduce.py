#!/usr/bin/env python3
# Copyright (c) d-SEAMS developers
"""Run every Reproduce the paper row, or the full rows when paper_full is on."""

import os
import shlex
import subprocess
import sys
from pathlib import Path

BEGIN = "# BEGIN PAPER_REPRODUCE"
END = "# END PAPER_REPRODUCE"


def fail(message):
    print(message, file=sys.stderr)
    raise SystemExit(1)


def entries(org):
    text = org.read_text(encoding="utf-8")
    start = text.find(BEGIN)
    stop = text.find(END)
    if start < 0 or stop < start:
        fail("reproduce.org is missing the paper catalog")
    rows = []
    for line in text[start:stop].splitlines():
        line = line.strip()
        if not line.startswith("|") or line.startswith("|-"):
            continue
        cells = [cell.strip() for cell in line.strip("|").split("|")]
        if cells[0] == "entry":
            continue
        if len(cells) != 8:
            fail("catalog row has {} cells: {}".format(len(cells), line))
        rows.append(
            {
                "entry": cells[0],
                "size": cells[1],
                "input": cells[2],
                "command": cells[3],
                "expect": cells[4],
                "regenerator": cells[5],
                "csv": cells[6],
                "csv_expect": cells[7],
            }
        )
    if len(rows) < 20:
        fail("catalog has {} rows".format(len(rows)))
    return rows


def check_regenerator(root, item):
    spec = item["regenerator"]
    if ":" not in spec:
        fail("{} regenerator needs file:snippet".format(item["entry"]))
    path, snippet = spec.split(":", 1)
    file = root / path
    if not file.is_file() or snippet not in file.read_text(encoding="utf-8"):
        fail("{} regenerator {} is missing".format(item["entry"], spec))


def check_csv(root, item):
    if not item["csv"]:
        return
    file = root / item["csv"]
    if not file.is_file():
        fail("{} csv is missing: {}".format(item["entry"], item["csv"]))
    if item["csv_expect"] and item["csv_expect"] not in file.read_text(encoding="utf-8"):
        fail("{} csv lacks {!r}".format(item["entry"], item["csv_expect"]))


def run(cmd, cwd, env):
    done = subprocess.run(cmd, cwd=str(cwd), env=env, capture_output=True, text=True)
    if done.returncode != 0:
        detail = (done.stdout or "") + (done.stderr or "")
        fail("exit {} for {}\n{}".format(done.returncode, " ".join(cmd), detail[-800:]))
    return done.stdout


def execute(root, seams, libdir, build, item, env):
    command = item["command"]
    out = Path(build) / "paper-reproduce" / item["entry"]
    out.mkdir(parents=True, exist_ok=True)
    if command.startswith("paper "):
        topic = command.split(None, 1)[1]
        script = root / "docs" / "tutorials" / "paper_science.py"
        return run([sys.executable, str(script), topic, seams, libdir, str(out)], root, env)
    if command.startswith("seams "):
        return run([seams, *shlex.split(command[len("seams "):])], root, env)
    fail("{} command is not paper or seams: {}".format(item["entry"], command))


def main():
    full = "--full" in sys.argv
    args = [item for item in sys.argv[1:] if item != "--full"]
    if len(args) != 4:
        fail("usage: check_paper_reproduce.py [--full] ROOT SEAMS LIBDIR BUILD")
    if full and os.environ.get("SEAMS_PAPER_FULL") != "1":
        print("paper_full is off")
        raise SystemExit(77)
    root, seams, libdir, build = args
    root = Path(root)
    env = os.environ.copy()
    env.setdefault("QT_QPA_PLATFORM", "offscreen")
    rows = entries(root / "docs" / "orgmode" / "paper" / "reproduce.org")
    wanted = "full" if full else "small"
    ran = 0
    for item in rows:
        check_regenerator(root, item)
        check_csv(root, item)
        if item["size"] != wanted:
            continue
        if item["size"] not in {"small", "full"}:
            fail("{} size is {}".format(item["entry"], item["size"]))
        text = execute(root, seams, libdir, build, item, env)
        missing = [part for part in item["expect"].split(" && ") if part not in text]
        if missing:
            fail("{} missing {}\n{}".format(item["entry"], missing, text[-500:]))
        ran += 1
        print(item["entry"], "ok")
    if ran < 1:
        fail("no {} rows".format(wanted))
    print("paper reproduce {} rows {}".format(wanted, ran))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
