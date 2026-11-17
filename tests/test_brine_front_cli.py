"""Regression check for brine_front multi-frame dumps.

usage: test_brine_front_cli.py BRINE_FRONT
"""

import subprocess
import sys


def main() -> int:
    brine_front = sys.argv[1]
    result = subprocess.run(
        [brine_front, "input/traj/nucleation_three.lammpstrj"],
        capture_output=True,
        text=True,
    )
    if result.returncode != 0:
        sys.stderr.write(
            "brine_front failed\n"
            + f"returncode={result.returncode}\n"
            + result.stdout
            + result.stderr
        )
        return 1

    stdout_lines = [line.strip() for line in result.stdout.splitlines() if line.strip()]
    frame_lines = [line for line in stdout_lines if line.startswith("frame ")]
    summary_lines = [line for line in stdout_lines if line.startswith("summary ")]

    if len(frame_lines) != 3:
        sys.stderr.write("wrong frame count\n" + result.stdout + result.stderr)
        return 1
    if frame_lines[0].split()[:2] != ["frame", "0"]:
        sys.stderr.write("first frame index\n" + frame_lines[0] + "\n")
        return 1
    if frame_lines[1].split()[:2] != ["frame", "1"]:
        sys.stderr.write("second frame index\n" + frame_lines[1] + "\n")
        return 1
    if frame_lines[2].split()[:2] != ["frame", "2"]:
        sys.stderr.write("third frame index\n" + frame_lines[2] + "\n")
        return 1

    expected_summary = (
        "summary frames 3 front_A_per_ns 0 front_stderr 0 thickness_slope 0 "
        "traps 0 releases 0 k_mean 0 c_liquid_last 0"
    )
    if summary_lines != [expected_summary]:
        sys.stderr.write("wrong summary\n" + result.stdout + result.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
