#!/usr/bin/env python3
import pathlib
import subprocess
import sys


def run(seams: str, trajectory: pathlib.Path) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        [
            seams,
            "--frame",
            "1",
            "--last",
            "3",
            "--type",
            "8",
            "--exclude-resname",
            "SOLX",
            "--exclude-atomname",
            "MW",
            "--type-map",
            "O_2=21,N_2=22,C_3=23,C_2=24,C_R=25,H_=26",
            "chill-plus",
            str(trajectory),
        ],
        check=False,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )


def main() -> int:
    if len(sys.argv) != 3:
        raise SystemExit("usage: test_confined_water_cli.py SEAMS FIXTURE-PREFIX")
    seams = sys.argv[1]
    prefix = pathlib.Path(sys.argv[2])
    gro = run(seams, prefix.with_suffix(".gro"))
    xtc = run(seams, prefix.with_suffix(".xtc"))
    if gro.returncode != 0:
        sys.stderr.write(gro.stderr)
        return gro.returncode
    if xtc.returncode != 0:
        sys.stderr.write(xtc.stderr)
        return xtc.returncode
    if gro.stdout.count("nop ") != 3:
        sys.stderr.write("GRO run did not produce three frames\n")
        return 1
    if xtc.stdout != gro.stdout:
        sys.stderr.write("XTC and GRO CHILL+ output differs\n")
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
