#!/usr/bin/env python3
# Copyright (c) d-SEAMS developers
"""ensure_cubic_trajectory() must not trust a pre-existing file on sight.

bench_h2h_rings.py's ensure_cubic_trajectory() used to return as soon as
input/traj/mW_cubic.lammpstrj existed, bypassing input/fetch.py's pinned
SHA-256 check for every pre-existing file. This test plants a corrupt
file at that path and checks that ensure_cubic_trajectory() still goes
through fetch.py's checksum verification and raises, instead of silently
using the corrupt file.
"""

from __future__ import annotations

import importlib.util
import io
import shutil
import sys
import tempfile
from pathlib import Path
from unittest import mock

ROOT = Path(__file__).resolve().parents[1]
CORRUPT_CONTENT = b"not the real trajectory"


class _FakeResponse:
    """Stands in for the object urllib.request.urlopen() returns."""

    def __init__(self, payload: bytes):
        self._stream = io.BytesIO(payload)

    def __enter__(self):
        return self._stream

    def __exit__(self, *exc_info):
        return False


def load_bench_h2h_rings(fake_root: Path):
    path = ROOT / "tests" / "bench_h2h_rings.py"
    spec = importlib.util.spec_from_file_location(
        "bench_h2h_rings_fetch_test", path
    )
    module = importlib.util.module_from_spec(spec)
    # dataclass field resolution looks the defining module up by name in
    # sys.modules, so it must be registered before exec_module runs.
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    module.ROOT = fake_root
    return module


def main() -> int:
    with tempfile.TemporaryDirectory(prefix="bench-h2h-fetch-test-") as work:
        fake_root = Path(work)
        (fake_root / "input").mkdir()
        shutil.copyfile(
            ROOT / "input" / "fetch.py", fake_root / "input" / "fetch.py"
        )
        traj_dir = fake_root / "input" / "traj"
        traj_dir.mkdir()
        corrupt = traj_dir / "mW_cubic.lammpstrj"
        corrupt.write_bytes(CORRUPT_CONTENT)

        bench = load_bench_h2h_rings(fake_root)

        # The mocked response also does not match the pinned digest, so a
        # checksum mismatch - not a network failure - is what must surface.
        with mock.patch(
            "urllib.request.urlopen",
            return_value=_FakeResponse(b"still not the real trajectory"),
        ):
            try:
                bench.ensure_cubic_trajectory()
            except RuntimeError as exc:
                message = str(exc)
                if "expected sha256" not in message:
                    print(f"wrong failure mode: {message}", file=sys.stderr)
                    return 1
                print(f"ensure_cubic_trajectory raised as expected: {message}")
            else:
                print(
                    "ensure_cubic_trajectory silently accepted a corrupt "
                    "pre-existing file instead of verifying it",
                    file=sys.stderr,
                )
                return 1

        if corrupt.read_bytes() != CORRUPT_CONTENT:
            print(
                "the corrupt file was replaced instead of left alone on "
                "a failed verification",
                file=sys.stderr,
            )
            return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
