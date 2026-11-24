#!/usr/bin/env python3
"""The versioned-fragment check accepts towncrier's continuation indent only."""

import tempfile
import unittest
from pathlib import Path
from unittest import mock

import check_changelog
from check_changelog import ROOT, check, continuation_bullet, missing_versioned


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
        self.assertEqual(len(versioned), 58)
        self.assertEqual(missing_versioned(changelog, fragments), {})


class TitleCheck(unittest.TestCase):
    def test_missing_bullet_still_reports_a_stale_title(self):
        """A missing bullet must not reuse the meson version or the name map.

        The org title is the stale 2.14.1 line. The fragment's second
        paragraph is not indented, so the bullet is missing. Both have to
        be reported, and the title error has to name the meson version.
        """
        body = (
            "The pixi full-suite now fails whenever a test is skipped.\n"
            "\n"
            "`-Dwith_lammps=enabled` now requires the headers.\n"
        )
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            release = root / "changelog.d" / "2.14.1"
            release.mkdir(parents=True)
            (release / "+lammps.changed.md").write_text(body, encoding="utf-8")
            (root / "CHANGELOG.md").write_text(
                "## [2.14.1] - 2026-11-16\n\n"
                "- The pixi full-suite now fails whenever a test is skipped.\n"
                "\n"
                "`-Dwith_lammps=enabled` now requires the headers.\n",
                encoding="utf-8",
            )
            page = root / "release-2.14.org"
            page.write_text(
                "#+TITLE: 2.14.1 release notes\n"
                ".. include:: release-2.14-notes.rst\n",
                encoding="utf-8",
            )
            with (
                mock.patch.object(check_changelog, "ROOT", root),
                mock.patch.object(check_changelog, "PAGE", page),
                mock.patch.object(
                    check_changelog, "FRAGMENTS", root / "changelog.d", create=True
                ),
                mock.patch.object(check_changelog, "TYPES", ("fixed",), create=True),
                mock.patch.object(check_changelog, "HEADINGS", {}, create=True),
            ):
                errors = check([], "", "2.14.2")
        text = "\n".join(errors)
        self.assertIn("missing the towncrier bullet", text)
        self.assertIn("+lammps.changed.md", text)
        self.assertIn(
            "title is '2.14.1 release notes', expected '2.14.2 release notes'",
            text,
        )
        self.assertNotIn("expected '2.14.1 release notes'", text)


if __name__ == "__main__":
    unittest.main()
