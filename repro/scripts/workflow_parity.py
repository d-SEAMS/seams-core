#!/usr/bin/env python3
"""Compare the stable d-SEAMS workflows across CLI, Python, and Lua."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import subprocess
import sys
from typing import Any


ICE_LABELS = (
    "cubic",
    "hexagonal",
    "water",
    "interfacial",
    "clathrate",
    "interClathrate",
    "reCubic",
    "reHex",
    "unclassified",
)
ARRAY_KEYS = {"rdf.r", "rdf.g", "density.centres", "density.rho"}
INTEGER_KEYS = (
    {
        "read.water_nop",
        "read.ice_nop",
        "chill.nop",
        "chill_plus.nop",
        "cages.nop",
        "cages.ih",
        "cages.ic",
        "cages.water",
        "hbonds.nop",
        "hbonds.count",
        "pairs.count",
        "pairs.n_cation",
        "pairs.n_anion",
        "domain.n",
        "domain.largest",
    }
    | {f"chill.{label}" for label in ICE_LABELS}
    | {f"chill_plus.{label}" for label in ICE_LABELS}
)


def _after(tokens: list[str], name: str, convert: type = float) -> Any:
    try:
        value = tokens[tokens.index(name) + 1]
    except (ValueError, IndexError) as exc:
        raise ValueError(f"missing {name!r} in {' '.join(tokens)!r}") from exc
    return convert(value)


def _ice_counts(command: str, text: str) -> dict[str, int]:
    prefix = command.replace("-", "_")
    tokens = text.split()
    values = {f"{prefix}.{label}": 0 for label in ICE_LABELS}
    values[f"{prefix}.nop"] = _after(tokens, "nop", int)
    for index in range(2, len(tokens), 2):
        label = tokens[index]
        if label in ICE_LABELS and index + 1 < len(tokens):
            values[f"{prefix}.{label}"] = int(tokens[index + 1])
    return values


def parse_cli_text(command: str, text: str) -> dict[str, Any]:
    """Normalize one CLI text payload into the shared result keys."""
    if command == "read":
        return {"read.nop": _after(text.split(), "nop", int)}
    if command in {"chill", "chill-plus"}:
        return _ice_counts(command, text)
    if command == "cages":
        tokens = text.split()
        return {
            "cages.nop": _after(tokens, "nop", int),
            "cages.ih": _after(tokens, "hexagonal", int),
            "cages.ic": _after(tokens, "cubic", int),
            "cages.water": _after(tokens, "water", int),
        }
    if command == "rdf":
        rows = [
            line.split()
            for line in text.splitlines()
            if line and not line.startswith("#")
        ]
        return {
            "rdf.r": [float(row[0]) for row in rows],
            "rdf.g": [float(row[1]) for row in rows],
        }
    if command == "cn":
        header = next(
            (line.split() for line in text.splitlines() if line.startswith("# types")),
            [],
        )
        return {"cn.value": _after(header, "cn", float)}
    if command == "hbonds":
        tokens = text.split()
        return {
            "hbonds.nop": _after(tokens, "nop", int),
            "hbonds.count": _after(tokens, "hbonds", int),
        }
    if command == "density-z":
        rows = [
            line.split()
            for line in text.splitlines()
            if line and not line.startswith("#")
        ]
        return {
            "density.centres": [float(row[0]) for row in rows],
            "density.rho": [float(row[1]) for row in rows],
        }
    if command == "pairs":
        tokens = text.split()
        return {
            "pairs.count": _after(tokens, "count", int),
            "pairs.n_cation": _after(tokens, "nCation", int),
            "pairs.n_anion": _after(tokens, "nAnion", int),
        }
    if command == "domains":
        tokens = text.split()
        return {
            "domain.n": _after(tokens, "n", int),
            "domain.largest": _after(tokens, "largest", int),
            "domain.percolation": _after(tokens, "P", float),
        }
    raise ValueError(f"unsupported CLI command {command!r}")


def _run(
    command: list[str], *, env: dict[str, str] | None = None
) -> subprocess.CompletedProcess[str]:
    process = subprocess.run(command, capture_output=True, text=True, env=env)
    if process.returncode != 0:
        rendered = " ".join(command)
        raise RuntimeError(
            f"command failed ({process.returncode}): {rendered}\n"
            f"stdout:\n{process.stdout}\nstderr:\n{process.stderr}"
        )
    return process


def _cli_command(
    binary: Path, command: str, path: Path, *options: str
) -> dict[str, Any]:
    env = os.environ.copy()
    env.update({"NO_COLOR": "1", "LC_ALL": "C"})
    process = _run(
        [
            str(binary),
            command,
            str(path),
            *options,
            "--format",
            "json",
            "--strict-input",
        ],
        env=env,
    )
    lines = [line for line in process.stdout.splitlines() if line.strip()]
    if len(lines) != 1:
        raise ValueError(f"{command} emitted {len(lines)} JSONL records")
    payload = json.loads(lines[0])
    if payload.get("schema") != "dseams.cli/v1" or payload.get("status") != 0:
        raise ValueError(f"invalid {command} envelope: {payload!r}")
    return parse_cli_text(command, payload["text"])


def collect_cli(binary: Path, water: Path, ice: Path, ions: Path) -> dict[str, Any]:
    values: dict[str, Any] = {}
    water_read = _cli_command(binary, "read", water, "--type", "2")
    ice_read = _cli_command(binary, "read", ice, "--type", "1")
    values["read.water_nop"] = water_read["read.nop"]
    values["read.ice_nop"] = ice_read["read.nop"]
    values.update(_cli_command(binary, "chill", ice, "--type", "1", "--cutoff", "3.5"))
    values.update(
        _cli_command(binary, "chill-plus", ice, "--type", "1", "--cutoff", "3.5")
    )
    values.update(
        _cli_command(
            binary,
            "cages",
            ice,
            "--type",
            "1",
            "--cutoff",
            "3.5",
            "--graph",
            "seeded",
            "-k",
            "4",
        )
    )
    values.update(
        _cli_command(
            binary, "rdf", water, "--types", "2,2", "--cutoff", "6", "--bins", "60"
        )
    )
    values.update(
        _cli_command(
            binary, "cn", water, "--types", "2,2", "--cutoff", "6", "--bins", "60"
        )
    )
    values.update(
        _cli_command(
            binary,
            "hbonds",
            water,
            "--type",
            "2",
            "--htype",
            "1",
            "--cutoff",
            "3.5",
        )
    )
    values.update(
        _cli_command(
            binary,
            "density-z",
            water,
            "--type",
            "2",
            "--axis",
            "z",
            "--bins",
            "20",
        )
    )
    values.update(
        _cli_command(
            binary,
            "pairs",
            ions,
            "--site",
            "1=cationHead,2=anion",
        )
    )
    values.update(
        _cli_command(
            binary,
            "domains",
            ions,
            "--site",
            "1=polar,2=apolar",
            "--subset",
            "polar",
            "--cutoff",
            "1.5",
        )
    )
    return values


def _cage_counts(hc: Any, ddc: Any) -> dict[str, int]:
    pairs = [(bool(h), bool(d)) for h, d in zip(hc, ddc)]
    return {
        "cages.nop": len(pairs),
        "cages.ih": sum(h for h, _ in pairs),
        "cages.ic": sum(d and not h for h, d in pairs),
        "cages.water": sum(not h and not d for h, d in pairs),
    }


def _hbond_count(network: Any) -> int:
    return sum(max(0, len(row) - 1) for row in network) // 2


def collect_python(water: Path, ice: Path, ions: Path) -> tuple[dict[str, Any], str]:
    import pydseams as ds

    ice_frame = ds.read(ice, atom_type=1, bonded="cutoff", cutoff=3.5)
    water_frame = ds.read(water, atom_type=2, bonded="hbond", cutoff=3.5)
    ions_frame = ds.read(ions, all_atoms=True, atom_type=1, bonded="cutoff")
    values: dict[str, Any] = {
        "read.water_nop": water_frame.n_atoms,
        "read.ice_nop": ice_frame.n_atoms,
    }
    chill = ice_frame.chill()
    values.update({f"chill.{label}": int(chill[label]) for label in ICE_LABELS})
    values["chill.nop"] = ice_frame.n_atoms
    chill_plus = ice_frame.chill_plus()
    values.update(
        {f"chill_plus.{label}": int(chill_plus[label]) for label in ICE_LABELS}
    )
    values["chill_plus.nop"] = ice_frame.n_atoms
    cages = ice_frame.cages(seeded=True, k=4, candidate_cutoff=5.0)
    values.update(_cage_counts(cages.hc, cages.ddc))
    rdf_r, rdf_g = water_frame.rdf(2, 2, cutoff=6.0, binwidth=0.1)
    values["rdf.r"] = [float(value) for value in rdf_r]
    values["rdf.g"] = [float(value) for value in rdf_g]
    values["cn.value"] = float(water_frame.cn(2, 2, cutoff=6.0, binwidth=0.1))
    values["hbonds.nop"] = water_frame.n_atoms
    values["hbonds.count"] = _hbond_count(water_frame.hbonds)
    density = water_frame.density(bins=20, axis="z", atom_type=2)
    values["density.centres"] = [float(value) for value in density.centres]
    values["density.rho"] = [float(value) for value in density.rho]
    pair_table = ds.yoda.parseSiteSpec("1=cationHead,2=anion")
    pairs = ions_frame.pairs(pair_table)
    values.update(
        {
            "pairs.count": pairs.count,
            "pairs.n_cation": pairs.n_cation,
            "pairs.n_anion": pairs.n_anion,
        }
    )
    domain_table = ds.yoda.parseSiteSpec("1=polar,2=apolar")
    domain = ions_frame.domain(domain_table, ds.yoda.Kind.polar, cutoff=1.5)
    values.update(
        {
            "domain.n": domain.n,
            "domain.largest": domain.largest,
            "domain.percolation": domain.percolation,
        }
    )
    return values, str(ds.__version__)


def collect_lua(
    lua: Path,
    script: Path,
    source: Path,
    build: Path,
    water: Path,
    ice: Path,
    ions: Path,
) -> tuple[dict[str, Any], str]:
    env = os.environ.copy()
    env["LUA_PATH"] = f"{source / 'lua' / '?.lua'};;"
    env["LUA_CPATH"] = f"{build / '?.so'};;"
    process = _run([str(lua), str(script), str(water), str(ice), str(ions)], env=env)
    values: dict[str, Any] = {}
    for line in process.stdout.splitlines():
        if not line.strip():
            continue
        try:
            key, raw = line.split("\t", 1)
        except ValueError as exc:
            raise ValueError(f"unexpected Lua output: {line!r}") from exc
        if key in ARRAY_KEYS:
            values[key] = [] if not raw else [float(value) for value in raw.split(",")]
        elif key in INTEGER_KEYS:
            values[key] = int(raw)
        else:
            values[key] = float(raw)
    version_process = _run([str(lua), "-v"])
    version = version_process.stdout.strip() or version_process.stderr.strip()
    return values, version


def _values_close(left: Any, right: Any) -> bool:
    if isinstance(left, list) or isinstance(right, list):
        if (
            not isinstance(left, list)
            or not isinstance(right, list)
            or len(left) != len(right)
        ):
            return False
        return all(_values_close(a, b) for a, b in zip(left, right))
    if isinstance(left, bool) or isinstance(right, bool):
        return left is right
    if isinstance(left, (int, float)) and isinstance(right, (int, float)):
        if isinstance(left, int) and isinstance(right, int):
            return left == right
        return math.isclose(float(left), float(right), rel_tol=2e-5, abs_tol=1e-8)
    return left == right


def compare_interfaces(interfaces: dict[str, dict[str, Any]]) -> list[dict[str, Any]]:
    """Compare every interface to the CLI and report each shared key."""
    if "cli" not in interfaces:
        raise ValueError("interfaces must include cli")
    names = sorted(interfaces)
    keys = sorted(set().union(*(values.keys() for values in interfaces.values())))
    checks: list[dict[str, Any]] = []
    for key in keys:
        present = {name: key in interfaces[name] for name in names}
        reference = interfaces["cli"].get(key)
        passed = all(present.values()) and all(
            _values_close(reference, interfaces[name][key])
            for name in names
            if name != "cli" and present[name]
        )
        checks.append({"key": key, "pass": passed, "present": present})
    return checks


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _version(command: list[str]) -> str:
    process = _run(command)
    return (process.stdout or process.stderr).strip()


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--seams", type=Path, required=True)
    parser.add_argument("--lua", type=Path, required=True)
    parser.add_argument("--lua-source", type=Path, required=True)
    parser.add_argument("--lua-build", type=Path, required=True)
    parser.add_argument("--water", type=Path, required=True)
    parser.add_argument("--ice", type=Path, required=True)
    parser.add_argument("--ions", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    paths = {
        name: getattr(args, name).resolve()
        for name in ("seams", "lua", "lua_source", "lua_build", "water", "ice", "ions")
    }
    for name, path in paths.items():
        if not path.exists():
            parser.error(f"--{name.replace('_', '-')} does not exist: {path}")

    cli = collect_cli(paths["seams"], paths["water"], paths["ice"], paths["ions"])
    python, python_version = collect_python(paths["water"], paths["ice"], paths["ions"])
    lua_script = Path(__file__).resolve().parent.parent / "lua" / "workflow_parity.lua"
    lua, lua_version = collect_lua(
        paths["lua"],
        lua_script,
        paths["lua_source"],
        paths["lua_build"],
        paths["water"],
        paths["ice"],
        paths["ions"],
    )
    interfaces = {"cli": cli, "python": python, "lua": lua}
    checks = compare_interfaces(interfaces)
    passed = all(check["pass"] for check in checks)
    report = {
        "schema": "dseams.workflow-parity/v1",
        "status": "pass" if passed else "fail",
        "inputs": {
            name: {"path": str(paths[name]), "sha256": _sha256(paths[name])}
            for name in ("water", "ice", "ions")
        },
        "versions": {
            "cli": _version([str(paths["seams"]), "--version"]),
            "python": python_version,
            "lua": lua_version,
        },
        "interfaces": interfaces,
        "checks": checks,
    }
    output = args.output.resolve()
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(
        json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(f"workflow parity: {report['status']} ({len(checks)} checks)")
    if not passed:
        for check in checks:
            if not check["pass"]:
                print(f"mismatch: {check['key']}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
