#!/usr/bin/env python3
# Copyright (c) d-SEAMS developers
"""Verify GitHub Actions workflow triggers and deploy guards for develop.

Every gating workflow (``ci_*.yml`` or any workflow with a ``pull_request``
trigger) must run on ``pull_request`` and ``push`` targeting both ``main`` and
``develop``. Deploy, publish, and release steps or workflows must never run on
``develop`` (neither push nor pull_request): either the workflow itself does
not trigger on ``develop`` (such as tag-only ``release.yml``), or the job/step
carries an explicit ``if:`` expression that evaluates to false for both
``push`` to ``refs/heads/develop`` and ``pull_request`` targeting ``develop``.

Usage:
  tests/check_workflow_triggers.py [WORKFLOWS_DIR]
  tests/check_workflow_triggers.py --self-test
"""

from __future__ import annotations

import re
import sys
import tempfile
from pathlib import Path
from typing import Any

sys.dont_write_bytecode = True

try:
    import yaml as _pyyaml
except ImportError:  # pragma: no cover - covered by stdlib parser self-test
    _pyyaml = None

ROOT = Path(__file__).resolve().parents[1]
WORKFLOWS_DIR = ROOT / ".github" / "workflows"

REPO_SLUG = "d-SEAMS/seams-core"

DEVELOP_CONTEXTS: dict[str, dict[str, str]] = {
    "push to develop": {
        "github.repository": REPO_SLUG,
        "github.event_name": "push",
        "github.ref": "refs/heads/develop",
        "github.ref_name": "develop",
        "github.ref_type": "branch",
        "github.base_ref": "",
    },
    "pull_request to develop": {
        "github.repository": REPO_SLUG,
        "github.event_name": "pull_request",
        "github.ref": "refs/pull/1/merge",
        "github.ref_name": "1/merge",
        "github.ref_type": "branch",
        "github.base_ref": "develop",
    },
}

DEPLOY_ACTION_RE = re.compile(
    r"(?:^|/)(?:actions-gh-pages|deploy-pages|gh-action-pypi-publish"
    r"|action-gh-release|release-action|[^/]*(?:deploy|publish|gh-pages|gh-release)[^/]*)$",
    re.IGNORECASE,
)

DEPLOY_RUN_PATTERNS = (
    re.compile(r"\bgh\s+release\s+(?:create|upload)\b"),
    re.compile(r"\bgit\s+push\b[^\n]*\bgh-pages\b"),
    re.compile(r"\btwine\s+upload\b"),
    re.compile(r"\bcargo\s+publish\b"),
    re.compile(r"\bnpm\s+publish\b"),
    re.compile(r"\banaconda\s+upload\b"),
    re.compile(r"\bconda\s+upload\b"),
    re.compile(r"\bdocker\s+push\b"),
)


# ---------------------------------------------------------------------------
# Minimal stdlib YAML parser for GitHub Actions workflow files (fallback when
# PyYAML is absent, e.g. on minimal or macOS test runners).
# ---------------------------------------------------------------------------


def _strip_comment(text: str) -> str:
    in_single = False
    in_double = False
    i = 0
    while i < len(text):
        ch = text[i]
        if in_single:
            if ch == "'":
                if i + 1 < len(text) and text[i + 1] == "'":
                    i += 2
                    continue
                in_single = False
        elif in_double:
            if ch == "\\" and i + 1 < len(text):
                i += 2
                continue
            if ch == '"':
                in_double = False
        else:
            if ch == "'":
                in_single = True
            elif ch == '"':
                in_double = True
            elif ch == "#" and (i == 0 or text[i - 1] in (" ", "\t")):
                return text[:i].rstrip()
        i += 1
    if in_single or in_double:
        raise ValueError(f"unterminated quoted string: {text!r}")
    return text.rstrip()


def _parse_scalar(raw: str) -> Any:
    s = raw.strip()
    if s == "" or s in ("null", "Null", "NULL", "~"):
        return None
    if s in ("true", "True", "TRUE"):
        return True
    if s in ("false", "False", "FALSE"):
        return False
    if s.startswith("'") and s.endswith("'") and len(s) >= 2:
        return s[1:-1].replace("''", "'")
    if s.startswith('"') and s.endswith('"') and len(s) >= 2:
        out: list[str] = []
        i = 1
        end = len(s) - 1
        while i < end:
            if s[i] == "\\" and i + 1 < end:
                nxt = s[i + 1]
                out.append({"n": "\n", "t": "\t", "r": "\r", "\\": "\\", '"': '"'}.get(nxt, nxt))
                i += 2
            else:
                out.append(s[i])
                i += 1
        return "".join(out)
    if s.startswith("["):
        if not s.endswith("]"):
            raise ValueError(f"unterminated flow sequence: {raw!r}")
        inner = s[1:-1].strip()
        if not inner:
            return []
        items: list[Any] = []
        buf: list[str] = []
        in_single = False
        in_double = False
        depth = 0
        i = 0
        while i < len(inner):
            ch = inner[i]
            if in_single:
                buf.append(ch)
                if ch == "'":
                    if i + 1 < len(inner) and inner[i + 1] == "'":
                        buf.append(inner[i + 1])
                        i += 2
                        continue
                    in_single = False
            elif in_double:
                buf.append(ch)
                if ch == "\\" and i + 1 < len(inner):
                    buf.append(inner[i + 1])
                    i += 2
                    continue
                if ch == '"':
                    in_double = False
            else:
                if ch == "'":
                    in_single = True
                    buf.append(ch)
                elif ch == '"':
                    in_double = True
                    buf.append(ch)
                elif ch in "[{":
                    depth += 1
                    buf.append(ch)
                elif ch in "]}":
                    depth -= 1
                    buf.append(ch)
                elif ch == "," and depth == 0:
                    items.append(_parse_scalar("".join(buf)))
                    buf = []
                else:
                    buf.append(ch)
            i += 1
        if depth != 0 or in_single or in_double:
            raise ValueError(f"malformed flow sequence: {raw!r}")
        if buf:
            items.append(_parse_scalar("".join(buf)))
        return items
    if s.startswith("{"):
        if not s.endswith("}"):
            raise ValueError(f"unterminated flow mapping: {raw!r}")
        inner = s[1:-1].strip()
        if not inner:
            return {}
        mapping: dict[str, Any] = {}
        for part in inner.split(","):
            k, sep, v = part.partition(":")
            if not sep:
                raise ValueError(f"malformed flow mapping entry: {part!r}")
            mapping[str(_parse_scalar(k))] = _parse_scalar(v)
        return mapping
    return s


def _split_key_value(line: str) -> tuple[str, str] | None:
    in_single = False
    in_double = False
    i = 0
    while i < len(line):
        ch = line[i]
        if in_single:
            if ch == "'":
                if i + 1 < len(line) and line[i + 1] == "'":
                    i += 2
                    continue
                in_single = False
        elif in_double:
            if ch == "\\" and i + 1 < len(line):
                i += 2
                continue
            if ch == '"':
                in_double = False
        else:
            if ch == "'":
                in_single = True
            elif ch == '"':
                in_double = True
            elif ch == ":" and (i + 1 == len(line) or line[i + 1] in (" ", "\t")):
                key_raw = line[:i].strip()
                val_raw = line[i + 1 :].strip()
                if not key_raw:
                    return None
                return str(_parse_scalar(key_raw)), val_raw
        i += 1
    return None


def _parse_yaml_stdlib(text: str) -> Any:
    raw_lines = text.splitlines()
    lines: list[tuple[int, str, str]] = []
    for lineno, raw in enumerate(raw_lines, start=1):
        stripped = raw.lstrip(" ")
        if stripped.startswith("\t"):
            raise ValueError(f"line {lineno}: tab indentation is not valid YAML")
        indent = len(raw) - len(stripped)
        lines.append((indent, stripped, raw))

    n = len(lines)

    def skip_blank_and_comments(idx: int) -> int:
        while idx < n:
            _, s, _ = lines[idx]
            if s == "" or s.startswith("#"):
                idx += 1
            else:
                break
        return idx

    def parse_block_scalar(idx: int, parent_indent: int, header: str) -> tuple[str, int]:
        fold = header.startswith(">")
        chomp = "clip"
        if len(header) > 1 and header[1] in ("-", "+"):
            chomp = "strip" if header[1] == "-" else "keep"
        block_indent: int | None = None
        collected: list[str] = []
        while idx < n:
            ind, s, raw = lines[idx]
            if s == "":
                collected.append("")
                idx += 1
                continue
            if block_indent is None:
                if ind <= parent_indent:
                    break
                block_indent = ind
            if ind < block_indent:
                break
            collected.append(raw[block_indent:])
            idx += 1
        while collected and collected[-1] == "" and chomp != "keep":
            collected.pop()
        if fold:
            value = " ".join(part.strip() for part in collected if part.strip())
        else:
            value = "\n".join(collected)
        if chomp == "clip" and value:
            value += "\n"
        return value, idx

    def parse_block(idx: int, min_indent: int) -> tuple[Any, int]:
        idx = skip_blank_and_comments(idx)
        if idx >= n or lines[idx][0] < min_indent:
            return None, idx
        ind, s, _ = lines[idx]
        if s == "-" or s.startswith("- "):
            return parse_seq(idx, ind)
        return parse_map(idx, ind)

    def parse_seq(idx: int, seq_indent: int) -> tuple[list[Any], int]:
        out: list[Any] = []
        while True:
            idx = skip_blank_and_comments(idx)
            if idx >= n:
                break
            ind, s, _ = lines[idx]
            if ind < seq_indent:
                break
            if ind != seq_indent or not (s == "-" or s.startswith("- ")):
                break
            rest = _strip_comment(s[1:].lstrip(" ")) if len(s) > 1 else ""
            if rest == "":
                val, idx = parse_block(idx + 1, seq_indent + 1)
                out.append(val)
                continue
            kv = _split_key_value(rest)
            if kv is not None:
                key, val_raw = kv
                item_indent = seq_indent + 2
                item: dict[str, Any] = {}
                idx += 1
                val, idx = parse_value(val_raw, idx, item_indent)
                item[key] = val
                more, idx = parse_map(idx, item_indent, existing=item)
                out.append(more)
            else:
                out.append(_parse_scalar(rest))
                idx += 1
        return out, idx

    def parse_value(val_raw: str, next_idx: int, key_indent: int) -> tuple[Any, int]:
        if val_raw in ("|", "|-", "|+", ">", ">-", ">+"):
            return parse_block_scalar(next_idx, key_indent, val_raw)
        if val_raw == "":
            peek = skip_blank_and_comments(next_idx)
            if peek < n and (
                lines[peek][0] > key_indent
                or (
                    lines[peek][0] == key_indent
                    and (lines[peek][1] == "-" or lines[peek][1].startswith("- "))
                )
            ):
                return parse_block(peek, lines[peek][0])
            return None, next_idx
        return _parse_scalar(val_raw), next_idx

    def parse_map(
        idx: int, map_indent: int, existing: dict[str, Any] | None = None
    ) -> tuple[dict[str, Any], int]:
        out: dict[str, Any] = {} if existing is None else existing
        while True:
            idx = skip_blank_and_comments(idx)
            if idx >= n:
                break
            ind, s, _ = lines[idx]
            if ind < map_indent:
                break
            if ind != map_indent:
                raise ValueError(f"line {idx + 1}: unexpected indentation {ind} (expected {map_indent})")
            if s == "-" or s.startswith("- "):
                break
            cleaned = _strip_comment(s)
            kv = _split_key_value(cleaned)
            if kv is None:
                raise ValueError(f"line {idx + 1}: expected mapping entry, got {s!r}")
            key, val_raw = kv
            val, idx = parse_value(val_raw, idx + 1, map_indent)
            out[key] = val
        return out, idx

    start = skip_blank_and_comments(0)
    if start >= n:
        return None
    result, end = parse_block(start, lines[start][0])
    end = skip_blank_and_comments(end)
    if end < n:
        raise ValueError(f"line {end + 1}: trailing content at indentation {lines[end][0]}")
    return result


def load_yaml(text: str, *, use_pyyaml: bool = True) -> Any:
    """Parse YAML via PyYAML when available, or the stdlib parser otherwise."""
    if use_pyyaml and _pyyaml is not None:
        return _pyyaml.safe_load(text)
    return _parse_yaml_stdlib(text)


# ---------------------------------------------------------------------------
# GitHub Actions trigger matching
# ---------------------------------------------------------------------------


def _workflow_on(doc: dict[str, Any]) -> dict[str, Any]:
    raw = doc.get("on")
    if raw is None and True in doc:
        # PyYAML implements YAML 1.1 where unquoted `on:` parses as True.
        raw = doc[True]
    if isinstance(raw, str):
        return {raw: None}
    if isinstance(raw, list):
        return {str(item): None for item in raw}
    if isinstance(raw, dict):
        return {str(k): v for k, v in raw.items()}
    raise ValueError("workflow is missing a valid 'on' mapping")


def _glob_to_regex(pattern: str) -> re.Pattern[str]:
    parts: list[str] = ["^"]
    i = 0
    while i < len(pattern):
        if pattern[i : i + 2] == "**":
            parts.append(".*")
            i += 2
        elif pattern[i] == "*":
            parts.append("[^/]*")
            i += 1
        elif pattern[i] == "?":
            parts.append("[^/]")
            i += 1
        else:
            parts.append(re.escape(pattern[i]))
            i += 1
    parts.append("$")
    return re.compile("".join(parts))


def _branch_matches_patterns(branch: str, patterns: list[Any], *, initial: bool) -> bool:
    matched = initial
    for raw in patterns:
        pat = str(raw)
        if pat.startswith("!"):
            if _glob_to_regex(pat[1:]).match(branch):
                matched = False
        else:
            if _glob_to_regex(pat).match(branch):
                matched = True
    return matched


def event_matches_branch(event_name: str, event_cfg: Any, branch: str) -> bool:
    """Return True if a GitHub Actions event filter matches ``branch``."""
    if event_cfg is None:
        return True
    if not isinstance(event_cfg, dict):
        return False
    has_branches = "branches" in event_cfg and event_cfg["branches"] is not None
    has_ignore = "branches-ignore" in event_cfg and event_cfg["branches-ignore"] is not None
    has_tags = ("tags" in event_cfg and event_cfg["tags"] is not None) or (
        "tags-ignore" in event_cfg and event_cfg["tags-ignore"] is not None
    )
    if has_branches:
        pats = event_cfg["branches"]
        if isinstance(pats, str):
            pats = [pats]
        return _branch_matches_patterns(branch, list(pats), initial=False)
    if has_ignore:
        pats = event_cfg["branches-ignore"]
        if isinstance(pats, str):
            pats = [pats]
        return not _branch_matches_patterns(branch, list(pats), initial=False)
    if event_name == "push" and has_tags:
        # A push trigger with only tags/tags-ignore never runs on branch pushes.
        return False
    return True


def workflow_triggers_on_branch(on_cfg: dict[str, Any], event: str, branch: str) -> bool:
    if event not in on_cfg:
        return False
    return event_matches_branch(event, on_cfg[event], branch)


# ---------------------------------------------------------------------------
# GitHub Actions `if:` expression evaluator (fail-closed three-valued logic)
# ---------------------------------------------------------------------------

_TOKEN_RE = re.compile(
    r"\s*(?:(&&|\|\||==|!=|[!(),])|'((?:[^']|'')*)'|([A-Za-z_][A-Za-z0-9_.-]*))"
)


def _tokenize_expr(expr: str) -> list[tuple[str, str]]:
    text = expr.strip()
    if text.startswith("${{") and text.endswith("}}"):
        text = text[3:-2].strip()
    tokens: list[tuple[str, str]] = []
    pos = 0
    while pos < len(text):
        match = _TOKEN_RE.match(text, pos)
        if match is None:
            raise ValueError(f"unsupported token in expression near {text[pos:]!r}")
        op, quoted, ident = match.groups()
        if op is not None:
            tokens.append(("OP", op))
        elif quoted is not None:
            tokens.append(("STR", quoted.replace("''", "'")))
        else:
            tokens.append(("IDENT", ident))
        pos = match.end()
    return tokens


def eval_if_expression(expr: Any, ctx: dict[str, str]) -> bool | None:
    """Evaluate a GitHub Actions ``if:`` expression against ``ctx``.

    Returns ``True``, ``False``, or ``None`` if the result cannot be statically
    determined from ``ctx``. Callers checking deploy guards must treat ``None``
    as unsafe (fail-closed).
    """
    if expr is None:
        return True
    if isinstance(expr, bool):
        return expr
    text = str(expr).strip()
    if not text:
        return True
    try:
        tokens = _tokenize_expr(text)
    except ValueError:
        return None
    if not tokens:
        return None

    idx = 0

    def peek() -> tuple[str, str] | None:
        return tokens[idx] if idx < len(tokens) else None

    def consume() -> tuple[str, str]:
        nonlocal idx
        tok = tokens[idx]
        idx += 1
        return tok

    def to_bool(val: Any) -> bool | None:
        if val is None:
            return None
        if isinstance(val, bool):
            return val
        if isinstance(val, str):
            return val != ""
        return bool(val)

    def parse_or() -> Any:
        left = parse_and()
        while peek() == ("OP", "||"):
            consume()
            right = parse_and()
            lb, rb = to_bool(left), to_bool(right)
            if lb is True or rb is True:
                left = True
            elif lb is False and rb is False:
                left = False
            else:
                left = None
        return left

    def parse_and() -> Any:
        left = parse_not()
        while peek() == ("OP", "&&"):
            consume()
            right = parse_not()
            lb, rb = to_bool(left), to_bool(right)
            if lb is False or rb is False:
                left = False
            elif lb is True and rb is True:
                left = True
            else:
                left = None
        return left

    def parse_not() -> Any:
        if peek() == ("OP", "!"):
            consume()
            val = to_bool(parse_not())
            return None if val is None else (not val)
        return parse_cmp()

    def parse_cmp() -> Any:
        left = parse_primary()
        tok = peek()
        if tok in (("OP", "=="), ("OP", "!=")):
            op = consume()[1]
            right = parse_primary()
            if left is None or right is None:
                return None
            if isinstance(left, str) and isinstance(right, str):
                eq = left.lower() == right.lower()
            else:
                eq = left == right
            return eq if op == "==" else (not eq)
        return left

    def parse_primary() -> Any:
        tok = peek()
        if tok is None:
            raise ValueError("unexpected end of expression")
        kind, val = consume()
        if kind == "OP" and val == "(":
            inner = parse_or()
            if peek() != ("OP", ")"):
                raise ValueError("missing closing ')'")
            consume()
            return inner
        if kind == "STR":
            return val
        if kind == "IDENT":
            lower = val.lower()
            if peek() == ("OP", "("):
                consume()
                args: list[Any] = []
                if peek() != ("OP", ")"):
                    args.append(parse_or())
                    while peek() == ("OP", ","):
                        consume()
                        args.append(parse_or())
                if peek() != ("OP", ")"):
                    raise ValueError("missing ')' in function call")
                consume()
                return eval_func(lower, args)
            if lower == "true":
                return True
            if lower == "false":
                return False
            return ctx.get(val)
        raise ValueError(f"unexpected token {val!r}")

    def eval_func(name: str, args: list[Any]) -> Any:
        if name == "always" and not args:
            return True
        if name == "success" and not args:
            return True
        if name in ("failure", "cancelled") and not args:
            return False
        if name == "startswith" and len(args) == 2:
            if isinstance(args[0], str) and isinstance(args[1], str):
                return args[0].lower().startswith(args[1].lower())
            return None
        if name == "endswith" and len(args) == 2:
            if isinstance(args[0], str) and isinstance(args[1], str):
                return args[0].lower().endswith(args[1].lower())
            return None
        if name == "contains" and len(args) == 2:
            if isinstance(args[0], str) and isinstance(args[1], str):
                return args[1].lower() in args[0].lower()
            return None
        return None

    try:
        result = to_bool(parse_or())
    except ValueError:
        return None
    if idx != len(tokens):
        return None
    return result


# ---------------------------------------------------------------------------
# Deploy / publish step detection and workflow verification
# ---------------------------------------------------------------------------


def deploy_step_reason(step: dict[str, Any]) -> str | None:
    """Return a short reason if ``step`` deploys, publishes, or cuts a release."""
    uses = step.get("uses")
    if isinstance(uses, str):
        action = uses.split("@", 1)[0].strip()
        if DEPLOY_ACTION_RE.search(action):
            return f"uses {action}"
    run = step.get("run")
    if isinstance(run, str):
        for pat in DEPLOY_RUN_PATTERNS:
            match = pat.search(run)
            if match is not None:
                return f"runs {match.group(0)!r}"
    return None


def _is_gating_workflow(path: Path, on_cfg: dict[str, Any]) -> bool:
    return path.name.startswith("ci_") or "pull_request" in on_cfg


def check_workflow_file(path: Path, *, use_pyyaml: bool = True) -> list[str]:
    problems: list[str] = []
    label = path.name
    try:
        doc = load_yaml(path.read_text(encoding="utf-8"), use_pyyaml=use_pyyaml)
    except Exception as exc:
        return [f"{label}: invalid YAML ({exc})"]
    if not isinstance(doc, dict):
        return [f"{label}: top-level YAML document must be a mapping"]

    try:
        on_cfg = _workflow_on(doc)
    except ValueError as exc:
        return [f"{label}: {exc}"]

    jobs = doc.get("jobs")
    if not isinstance(jobs, dict) or not jobs:
        return [f"{label}: missing non-empty 'jobs' mapping"]

    pr_develop = workflow_triggers_on_branch(on_cfg, "pull_request", "develop")
    pr_main = workflow_triggers_on_branch(on_cfg, "pull_request", "main")
    push_develop = workflow_triggers_on_branch(on_cfg, "push", "develop")
    push_main = workflow_triggers_on_branch(on_cfg, "push", "main")

    if _is_gating_workflow(path, on_cfg):
        if not pr_develop:
            problems.append(f"{label}: pull_request trigger does not include 'develop'")
        if not pr_main:
            problems.append(f"{label}: pull_request trigger does not include 'main'")
        if not push_develop:
            problems.append(f"{label}: push trigger does not include 'develop'")
        if not push_main:
            problems.append(f"{label}: push trigger does not include 'main'")

    active_develop_events: list[tuple[str, dict[str, str]]] = []
    if push_develop:
        active_develop_events.append(("push to develop", DEVELOP_CONTEXTS["push to develop"]))
    if pr_develop:
        active_develop_events.append(
            ("pull_request to develop", DEVELOP_CONTEXTS["pull_request to develop"])
        )

    for job_id, job in jobs.items():
        if not isinstance(job, dict):
            problems.append(f"{label}: job {job_id!r} is not a mapping")
            continue
        job_if = job.get("if")
        job_uses = job.get("uses")
        if isinstance(job_uses, str):
            action = job_uses.split("@", 1)[0].strip()
            if DEPLOY_ACTION_RE.search(action):
                for ctx_name, ctx in active_develop_events:
                    if eval_if_expression(job_if, ctx) is not False:
                        problems.append(
                            f"{label}: deploy job {job_id!r} (uses {action}) "
                            f"is not blocked on {ctx_name}"
                        )
        steps = job.get("steps") or []
        if not isinstance(steps, list):
            continue
        for step_index, step in enumerate(steps, start=1):
            if not isinstance(step, dict):
                continue
            reason = deploy_step_reason(step)
            if reason is None:
                continue
            step_name = step.get("name") or f"step #{step_index}"
            step_if = step.get("if")
            for ctx_name, ctx in active_develop_events:
                job_val = eval_if_expression(job_if, ctx)
                step_val = eval_if_expression(step_if, ctx)
                blocked = (job_val is False) or (step_val is False)
                if not blocked:
                    problems.append(
                        f"{label}: deploy step {job_id}/{step_name!r} ({reason}) "
                        f"is not blocked on {ctx_name}"
                    )

    return problems


def check_workflows_dir(workflows_dir: Path, *, use_pyyaml: bool = True) -> list[str]:
    if not workflows_dir.is_dir():
        return [f"{workflows_dir}: workflows directory does not exist"]
    files = sorted(
        p for p in workflows_dir.iterdir() if p.is_file() and p.suffix in (".yml", ".yaml")
    )
    if not files:
        return [f"{workflows_dir}: no workflow YAML files found"]
    problems: list[str] = []
    gating_count = 0
    for path in files:
        problems.extend(check_workflow_file(path, use_pyyaml=use_pyyaml))
        if path.name.startswith("ci_"):
            gating_count += 1
    if gating_count == 0:
        problems.append(f"{workflows_dir}: no gating ci_*.yml workflows found")
    return problems


def _extract_summary(doc: dict[str, Any]) -> dict[str, Any]:
    on_cfg = _workflow_on(doc)
    jobs_summary: dict[str, Any] = {}
    for job_id, job in (doc.get("jobs") or {}).items():
        steps_summary = []
        for step in job.get("steps") or []:
            if isinstance(step, dict):
                steps_summary.append(
                    (step.get("name"), step.get("if"), deploy_step_reason(step))
                )
        jobs_summary[job_id] = (job.get("if"), steps_summary)
    return {
        "triggers": {
            (ev, br): workflow_triggers_on_branch(on_cfg, ev, br)
            for ev in ("push", "pull_request")
            for br in ("main", "develop")
        },
        "jobs": jobs_summary,
    }


def self_test() -> list[str]:
    problems: list[str] = []

    # 1. The real repository workflows pass under both the active loader and
    #    the zero-dependency stdlib loader.
    real_problems = check_workflows_dir(WORKFLOWS_DIR, use_pyyaml=True)
    if real_problems:
        problems.extend(f"repo check failed: {p}" for p in real_problems)
    stdlib_problems = check_workflows_dir(WORKFLOWS_DIR, use_pyyaml=False)
    if stdlib_problems:
        problems.extend(f"repo stdlib check failed: {p}" for p in stdlib_problems)

    # 2. When PyYAML is installed, confirm stdlib parser agrees with PyYAML on
    #    every workflow file in .github/workflows.
    if _pyyaml is not None:
        for path in sorted(WORKFLOWS_DIR.glob("*.y*ml")):
            text = path.read_text(encoding="utf-8")
            py_sum = _extract_summary(load_yaml(text, use_pyyaml=True))
            std_sum = _extract_summary(load_yaml(text, use_pyyaml=False))
            if py_sum != std_sum:
                problems.append(f"{path.name}: stdlib parser disagreed with PyYAML")

    # 3. Synthetic positive and negative cases (tested with both parsers when
    #    PyYAML is available).
    cases: list[tuple[str, str, str, str | None]] = [
        (
            "gating workflow missing develop on pull_request",
            "ci_sample.yml",
            """
name: Sample
on:
  push:
    branches: [main, develop]
  pull_request:
    branches: [main]
jobs:
  test:
    runs-on: ubuntu-latest
    steps:
      - run: echo ok
""",
            "pull_request trigger does not include 'develop'",
        ),
        (
            "gating workflow missing develop on push",
            "ci_sample.yml",
            """
name: Sample
on:
  push:
    branches: [main]
  pull_request:
    branches: [main, develop]
jobs:
  test:
    runs-on: ubuntu-latest
    steps:
      - run: echo ok
""",
            "push trigger does not include 'develop'",
        ),
        (
            "gating workflow negating develop",
            "ci_sample.yml",
            """
name: Sample
on:
  push:
    branches: ['**', '!develop']
  pull_request:
    branches-ignore: [develop]
jobs:
  test:
    runs-on: ubuntu-latest
    steps:
      - run: echo ok
""",
            "pull_request trigger does not include 'develop'",
        ),
        (
            "ungated gh-pages deploy on develop",
            "ci_docs.yml",
            """
name: Docs
on:
  push:
    branches: [main, develop]
  pull_request:
    branches: [main, develop]
jobs:
  build-deploy:
    runs-on: ubuntu-latest
    steps:
      - name: Deploy to GitHub Pages
        uses: peaceiris/actions-gh-pages@84c30a85c19949d7eee79c4ff27748b70285e453
""",
            "not blocked on push to develop",
        ),
        (
            "deploy gated only by event_name == push (still fires on develop push)",
            "ci_docs.yml",
            """
name: Docs
on:
  push:
    branches: [main, develop]
  pull_request:
    branches: [main, develop]
jobs:
  build-deploy:
    runs-on: ubuntu-latest
    steps:
      - name: Deploy to GitHub Pages
        if: github.event_name == 'push'
        uses: peaceiris/actions-gh-pages@84c30a85c19949d7eee79c4ff27748b70285e453
""",
            "not blocked on push to develop",
        ),
        (
            "deploy explicitly targeting develop",
            "ci_docs.yml",
            """
name: Docs
on:
  push:
    branches: [main, develop]
  pull_request:
    branches: [main, develop]
jobs:
  build-deploy:
    runs-on: ubuntu-latest
    steps:
      - name: Deploy to GitHub Pages
        if: github.ref == 'refs/heads/develop'
        uses: peaceiris/actions-gh-pages@84c30a85c19949d7eee79c4ff27748b70285e453
""",
            "not blocked on push to develop",
        ),
        (
            "release workflow accidentally triggered on develop",
            "release.yml",
            """
name: Release
on:
  push:
    branches: [develop]
    tags: ["v*"]
jobs:
  publish:
    runs-on: ubuntu-latest
    steps:
      - name: Create release
        run: gh release create v1.0.0
""",
            "not blocked on push to develop",
        ),
        (
            "mixed docs workflow with main-only deploy guard passes",
            "ci_docs.yml",
            """
name: Docs
on:
  push:
    branches: [main, develop]
    tags: ["v*"]
  pull_request:
    branches: [main, develop]
jobs:
  build-deploy:
    runs-on: ubuntu-latest
    steps:
      - name: Deploy to GitHub Pages
        if: github.repository == 'd-SEAMS/seams-core' && github.event_name == 'push' && github.ref == 'refs/heads/main'
        uses: peaceiris/actions-gh-pages@84c30a85c19949d7eee79c4ff27748b70285e453
""",
            None,
        ),
        (
            "tag-only release workflow passes",
            "release.yml",
            """
name: Release
on:
  push:
    tags: ["v*"]
jobs:
  publish:
    runs-on: ubuntu-latest
    steps:
      - name: Create release
        run: gh release create v1.0.0
""",
            None,
        ),
        (
            "invalid YAML fails",
            "ci_broken.yml",
            "on:\n\tpush:\n\t\tbranches: [main, develop]\n",
            "invalid YAML",
        ),
    ]

    modes = [True, False] if _pyyaml is not None else [False]
    with tempfile.TemporaryDirectory(prefix="seams-wf-selftest-") as tmp:
        tmp_path = Path(tmp)
        for idx, (label, filename, body, want_err) in enumerate(cases):
            case_dir = tmp_path / f"case-{idx}"
            case_dir.mkdir()
            wf_path = case_dir / filename
            wf_path.write_text(body.lstrip("\n"), encoding="utf-8")
            for use_pyyaml in modes:
                mode_name = "pyyaml" if use_pyyaml else "stdlib"
                got = check_workflow_file(wf_path, use_pyyaml=use_pyyaml)
                if want_err is None:
                    if got:
                        problems.append(
                            f"{label} ({mode_name}): expected pass, got {'; '.join(got)}"
                        )
                else:
                    if not any(want_err in item for item in got):
                        problems.append(
                            f"{label} ({mode_name}): expected {want_err!r}, got {'; '.join(got) or 'PASS'}"
                        )

    return problems


def main(argv: list[str]) -> int:
    if argv[1:] == ["--self-test"]:
        problems = self_test()
        ok_msg = "workflow trigger check self-test ok"
    elif len(argv) <= 2:
        target = Path(argv[1]) if len(argv) == 2 else WORKFLOWS_DIR
        problems = check_workflows_dir(target)
        count = len(
            [p for p in target.iterdir() if p.is_file() and p.suffix in (".yml", ".yaml")]
        ) if target.is_dir() else 0
        ok_msg = f"workflow trigger check ok ({count} workflows)"
    else:
        print("usage: tests/check_workflow_triggers.py [--self-test | WORKFLOWS_DIR]", file=sys.stderr)
        return 2

    if problems:
        print("\n".join(problems), file=sys.stderr)
        return 1
    print(ok_msg)
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
