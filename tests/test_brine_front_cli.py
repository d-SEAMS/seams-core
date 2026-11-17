"""Regression check for brine_front multi-frame dumps.

usage: test_brine_front_cli.py BRINE_FRONT
"""

import subprocess
import sys


def _parse_fields(line: str) -> dict[str, float]:
    tokens = line.split()
    fields: dict[str, float] = {}
    for index in range(0, len(tokens), 2):
        key = tokens[index]
        value = tokens[index + 1]
        fields[key] = float(value)
    return fields


def _close(value: float, target: float, tol: float) -> bool:
    return abs(value - target) <= tol


def _linear_fit(x: list[float], y: list[float]) -> tuple[float, float]:
    n = len(x)
    x_mean = sum(x) / n
    y_mean = sum(y) / n
    sxx = sum((v - x_mean) ** 2 for v in x)
    sxy = sum((x[i] - x_mean) * (y[i] - y_mean) for i in range(n))
    slope = sxy / sxx if sxx > 0 else 0.0
    residual = sum((y[i] - (y_mean + slope * (x[i] - x_mean))) ** 2 for i in range(n))
    slope_stderr = (residual / (n - 2) / sxx) ** 0.5 if n > 2 and sxx > 0 else 0.0
    return slope, slope_stderr


def main() -> int:
    brine_front = sys.argv[1]
    result = subprocess.run(
        [brine_front, "input/traj/mW_cubic.lammpstrj"],
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

    if "You entered a frame that doesn't exist." in result.stdout or "You entered a frame that doesn't exist." in result.stderr:
        sys.stderr.write("unexpected out-of-range probe\n" + result.stdout + result.stderr)
        return 1
    if len(frame_lines) != 11:
        sys.stderr.write("wrong frame count\n" + result.stdout + result.stderr)
        return 1
    if frame_lines[0].split()[:2] != ["frame", "0"] or frame_lines[-1].split()[:2] != ["frame", "10"]:
        sys.stderr.write("frame indices\n" + frame_lines[0] + "\n" + frame_lines[-1] + "\n")
        return 1
    if len(summary_lines) != 1:
        sys.stderr.write("missing summary\n" + result.stdout + result.stderr)
        return 1

    # The first mW frame is a fully cubic slab: layer ring/cubic counts are
    # validated independently in the ring and interface-order test suites.
    frame0 = _parse_fields(frame_lines[0])
    if not _close(frame0["thickness"], 50.96, 0.1):
        sys.stderr.write("frame0 thickness\n" + frame_lines[0] + "\n")
        return 1
    if int(round(frame0["cubic"])) != 512 or int(round(frame0["rings6"])) != 256:
        sys.stderr.write("frame0 ring/cubic counts\n" + frame_lines[0] + "\n")
        return 1
    if int(round(frame0["hex"])) != 0:
        sys.stderr.write("frame0 hex count\n" + frame_lines[0] + "\n")
        return 1

    summary = _parse_fields(summary_lines[0].replace("summary ", "", 1))
    if int(round(summary["frames"])) != 11:
        sys.stderr.write("summary frame count\n" + summary_lines[0] + "\n")
        return 1
    times = []
    thickness = []
    for row in frame_lines:
        fields = _parse_fields(row)
        times.append(fields["t_ns"])
        thickness.append(fields["thickness"])
    slope, slope_stderr = _linear_fit(times, thickness)
    if not _close(summary["thickness_slope"], slope, 1e-6):
        sys.stderr.write("thickness slope mismatch\n" + summary_lines[0] + "\n")
        return 1
    if not _close(summary["front_A_per_ns"], slope / 2.0, 1e-6):
        sys.stderr.write("front velocity mismatch\n" + summary_lines[0] + "\n")
        return 1
    if not _close(summary["front_stderr"], slope_stderr / 2.0, 1e-6):
        sys.stderr.write("front stderr mismatch\n" + summary_lines[0] + "\n")
        return 1
    if summary["front_A_per_ns"] == 0.0:
        sys.stderr.write("zero front velocity\n" + summary_lines[0] + "\n")
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
