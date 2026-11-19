#!/usr/bin/env python3
"""Tests for reader-only text patterns."""

import unittest

from check_reader_text import PROCESS_PATTERNS, README_PATTERNS, source_pages
from registry_claims import ClaimParseError, claims


def labels(text, patterns):
    return {label for label, pattern in patterns if pattern.search(text)}


class ReaderTextPatterns(unittest.TestCase):
    def test_book_rejects_meson_test(self):
        self.assertIn(
            "documentation test command",
            labels("meson test -C bbdir", PROCESS_PATTERNS),
        )

    def test_book_rejects_docs_test_name(self):
        self.assertIn(
            "documentation test command",
            labels("docs_reader_text", PROCESS_PATTERNS),
        )

    def test_readme_allows_build_commands(self):
        self.assertNotIn(
            "documentation test command",
            labels("meson test -C bbdir", README_PATTERNS),
        )

    def test_rejects_registry_publication_status(self):
        self.assertIn(
            "registry publication status",
            labels("Registry publication may follow the release.", PROCESS_PATTERNS),
        )

    def test_rejects_luarocks_dseams(self):
        for text in ("LuaRocks ~dseams~", "LuaRocks package ``dseams``"):
            with self.subTest(text=text):
                self.assertIn(
                    "unpublished LuaRocks package",
                    labels(text, PROCESS_PATTERNS),
                )

    def test_rejects_fpm_registry_dseams(self):
        self.assertIn(
            "unpublished fpm package",
            labels("fpm package ``dseams``", PROCESS_PATTERNS),
        )
        self.assertNotIn(
            "unpublished fpm package",
            labels("fpm Git dependency ``dseams``", PROCESS_PATTERNS),
        )

    def test_scans_clean_language_index(self):
        pages = {rel: path for path, rel in source_pages()}
        self.assertIn("lang/index.rst", pages)
        text = pages["lang/index.rst"].read_text(encoding="utf-8")
        self.assertFalse(labels(text, PROCESS_PATTERNS))

    def test_rejects_build_status_badges(self):
        for text in (
            "![Build Status](badge.svg)",
            "![built with nix](badge.svg)",
            "Check our build status here",
        ):
            with self.subTest(text=text):
                self.assertIn(
                    "build status badge",
                    labels(text, README_PATTERNS),
                )

    def test_claim_parser_counts_expected_commands(self):
        line = (
            "pip install dseamslib==2.14.1; "
            "cargo add seams-core@2.14.1; "
            "conda install -c conda-forge seams-core=2.14.1"
        )
        found = list(claims(line))
        self.assertEqual(len(found), 3)

    def test_claim_parser_rejects_malformed_conda_line(self):
        with self.assertRaises(ClaimParseError):
            list(claims('conda install "seams-core=2.14.1'))


if __name__ == "__main__":
    unittest.main()
