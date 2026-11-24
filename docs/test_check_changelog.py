#!/usr/bin/env python3
"""The versioned-fragment check accepts towncrier's continuation indent only."""

import tempfile
import unittest
from pathlib import Path

from check_changelog import ROOT, continuation_bullet, missing_versioned


class ContinuationBullet(unittest.TestCase):
    def test_two_space_continuation_matches(self):
        body = "The first line wraps\nonto the next line.\n"
        changelog = (
            "## [2.14.1] - 2026-11-16\n\n"
            + continuation_bullet(body)
            + "\n\n## [2.14.0] - 2026-11-15\n"
        )
        with tempfile.TemporaryDirectory() as tmp:
            fragment = Path(tmp) / "2.14.1"
            fragment.mkdir()
            (fragment / "+wrapped.fixed.md").write_text(body, encoding="utf-8")
            self.assertEqual(missing_versioned(changelog, Path(tmp)), {})

    def test_blank_line_paragraph_needs_the_indent(self):
        body = (
            "The pixi full-suite now fails whenever a test is skipped.\n"
            "\n"
            "`-Dwith_lammps=enabled` now requires the headers.\n"
        )
        indented = (
            "## [2.14.1] - 2026-11-16\n\n"
            + continuation_bullet(body)
            + "\n"
        )
        unindented = (
            "## [2.14.1] - 2026-11-16\n\n"
            "- The pixi full-suite now fails whenever a test is skipped.\n"
            "\n"
            "`-Dwith_lammps=enabled` now requires the headers.\n"
        )
        with tempfile.TemporaryDirectory() as tmp:
            fragment = Path(tmp) / "2.14.1"
            fragment.mkdir()
            (fragment / "+lammps.changed.md").write_text(body, encoding="utf-8")
            self.assertEqual(missing_versioned(indented, Path(tmp)), {})
            self.assertEqual(
                missing_versioned(unindented, Path(tmp)),
                {"2.14.1": ["+lammps.changed.md"]},
            )

    def test_internal_spaces_are_not_folded(self):
        body = "Alpha ends. Beta starts.\n"
        doubled = (
            "## [2.14.1] - 2026-11-16\n\n"
            "- Alpha ends.  Beta starts.\n"
        )
        with tempfile.TemporaryDirectory() as tmp:
            fragment = Path(tmp) / "2.14.1"
            fragment.mkdir()
            (fragment / "+spaces.fixed.md").write_text(body, encoding="utf-8")
            self.assertEqual(
                missing_versioned(doubled, Path(tmp)),
                {"2.14.1": ["+spaces.fixed.md"]},
            )

    def test_released_fragments_match_the_changelog(self):
        changelog = (ROOT / "CHANGELOG.md").read_text(encoding="utf-8")
        fragments = ROOT / "changelog.d"
        versioned = [
            path
            for directory in fragments.iterdir()
            if directory.is_dir()
            for path in directory.glob("*.md")
        ]
        self.assertEqual(len(versioned), 45)
        self.assertEqual(missing_versioned(changelog, fragments), {})


if __name__ == "__main__":
    unittest.main()
