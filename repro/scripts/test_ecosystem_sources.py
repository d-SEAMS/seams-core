"""Unit tests for the locked ecosystem source preparation."""

from __future__ import annotations

import importlib.util
import json
from pathlib import Path
import tempfile
import unittest


SCRIPT = Path(__file__).with_name("ecosystem_sources.py")
SPEC = importlib.util.spec_from_file_location("ecosystem_sources", SCRIPT)
if SPEC is None or SPEC.loader is None:
    raise RuntimeError(f"cannot load {SCRIPT}")
ecosystem_sources = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(ecosystem_sources)


def write_lock(path: Path, revision: str = "a" * 40) -> None:
    path.write_text(
        json.dumps(
            {
                "schema": "dseams.ecosystem-lock/v1",
                "components": {
                    "linkcell": {
                        "directory": "linkcell",
                        "repository": "https://example.test/linkcell.git",
                        "revision": revision,
                    },
                    "pydseamslib": {
                        "directory": "PydSEAMSlib",
                        "repository": "https://example.test/PydSEAMSlib.git",
                        "revision": "b" * 40,
                    },
                    "yodastruct": {
                        "directory": "yodaStruct",
                        "repository": "https://example.test/yodaStruct.git",
                        "revision": "c" * 40,
                    },
                },
            }
        ),
        encoding="utf-8",
    )


class LockTests(unittest.TestCase):
    def test_load_lock_accepts_full_revisions(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            lock_path = Path(directory) / "lock.json"
            write_lock(lock_path)
            lock = ecosystem_sources.load_lock(lock_path)
            self.assertEqual(lock["components"]["linkcell"]["revision"], "a" * 40)

    def test_load_lock_rejects_symbolic_revision(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            lock_path = Path(directory) / "lock.json"
            write_lock(lock_path, "main")
            with self.assertRaisesRegex(ValueError, "40-character"):
                ecosystem_sources.load_lock(lock_path)


class WireTests(unittest.TestCase):
    def test_wire_links_every_binding_to_one_core(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            lock_path = root / "lock.json"
            write_lock(lock_path)
            sources = root / "sources"
            core = root / "seams-core"
            for name in ("linkcell", "PydSEAMSlib", "yodaStruct"):
                (sources / name / "subprojects").mkdir(parents=True)
            (core / "subprojects").mkdir(parents=True)

            ecosystem_sources.wire(
                ecosystem_sources.load_lock(lock_path), sources, core
            )

            self.assertEqual(
                (core / "subprojects" / "linkcell").resolve(),
                (sources / "linkcell").resolve(),
            )
            for name in ("PydSEAMSlib", "yodaStruct"):
                self.assertEqual(
                    (sources / name / "subprojects" / "seams-core").resolve(),
                    core.resolve(),
                )

    def test_wire_refuses_to_replace_a_real_path(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            lock_path = root / "lock.json"
            write_lock(lock_path)
            sources = root / "sources"
            core = root / "seams-core"
            for name in ("linkcell", "PydSEAMSlib", "yodaStruct"):
                (sources / name / "subprojects").mkdir(parents=True)
            occupied = core / "subprojects" / "linkcell"
            occupied.mkdir(parents=True)

            with self.assertRaisesRegex(FileExistsError, "refusing to replace"):
                ecosystem_sources.wire(
                    ecosystem_sources.load_lock(lock_path), sources, core
                )


if __name__ == "__main__":
    unittest.main()
