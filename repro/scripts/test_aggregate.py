#!/usr/bin/env python3
"""Regression tests for the paper-manifest aggregator."""

from __future__ import annotations

import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest


SCRIPT = Path(__file__).with_name("aggregate.py")


class AggregateTests(unittest.TestCase):
    def test_source_and_parity_manifests_are_embedded(self):
        source_manifest = {
            "schema": "dseams.source-manifest/v1",
            "components": {"seams-core": {"revision": "a" * 40}},
        }
        parity_manifest = {
            "schema": "dseams.workflow-parity/v1",
            "status": "pass",
            "checks": {"read.water_nop": {"status": "pass"}},
        }
        with tempfile.TemporaryDirectory() as tmp:
            results = Path(tmp)
            (results / "source-manifest.json").write_text(
                json.dumps(source_manifest), encoding="utf-8"
            )
            (results / "workflow-parity.json").write_text(
                json.dumps(parity_manifest), encoding="utf-8"
            )
            process = subprocess.run(
                [sys.executable, str(SCRIPT), str(results)],
                check=True,
                capture_output=True,
                text=True,
            )

        paper_manifest = json.loads(process.stdout)
        self.assertEqual(paper_manifest["source_manifest"], source_manifest)
        self.assertEqual(paper_manifest["workflow_parity"], parity_manifest)


if __name__ == "__main__":
    unittest.main()
