#!/usr/bin/env python3
# Copyright (c) d-SEAMS developers
"""Read a vhs tape down to the commands a headless run should execute."""

import shlex


def _unquote(raw):
    parts = shlex.split(raw)
    if parts and parts[0] == "Type":
        return " ".join(parts[1:])
    return raw.strip()


def parse_tape(text):
    """Return commands, expect fragments, test names, lead, duration, still dump."""
    commands = []
    expects = []
    tests = []
    lead = None
    duration = None
    still = None
    hidden = False
    pending = None
    for raw in text.splitlines():
        line = raw.strip()
        if not line:
            continue
        if line.startswith("# expect "):
            expects.append(line[len("# expect ") :].strip())
            continue
        if line.startswith("# test "):
            tests.append(line[len("# test ") :].strip())
            continue
        if line.startswith("# lead "):
            lead = float(line.split()[2])
            continue
        if line.startswith("# duration "):
            duration = float(line.split()[2])
            continue
        if line.startswith("# still "):
            still = line[len("# still ") :].strip()
            continue
        if line == "Hide":
            hidden = True
            continue
        if line == "Show":
            hidden = False
            continue
        if line.startswith("Type "):
            pending = _unquote(line)
            continue
        if line == "Enter" and pending is not None:
            if not hidden and not pending.startswith("PS1="):
                commands.append(pending)
            pending = None
            continue
    if lead is None or duration is None:
        raise ValueError("tape is missing lead or duration")
    return {
        "commands": commands,
        "expects": expects,
        "tests": tests,
        "lead": lead,
        "duration": duration,
        "still": still,
    }


def quote_type(command):
    return "Type {}".format(shlex.quote(command))
