#!/usr/bin/env python3
"""Tests for reader-only text patterns."""

import contextlib
import io
import tempfile
import unittest
from pathlib import Path

from check_reader_text import (
    PROCESS_PATTERNS,
    README_PATTERNS,
    main,
    normalised_text,
    source_pages,
)
from registry_claims import ClaimParseError, claims


def labels(text, patterns):
    text = normalised_text(text)
    return {label for label, pattern in patterns if pattern.search(text)}


class ReaderTextPatterns(unittest.TestCase):
    def test_normalised_text_keeps_newlines(self):
        """Horizontal whitespace collapses. A newline stays a newline.

        Folding ``\\s+`` would turn the wrapped command into one line.
        """
        self.assertEqual(
            normalised_text("meson \t test\nwrapped"),
            "meson test\nwrapped",
        )
        self.assertEqual(normalised_text("alpha\n\nbeta"), "alpha\n\nbeta")

    def test_book_rejects_meson_test(self):
        self.assertIn(
            "documentation test command",
            labels("meson test -C bbdir", PROCESS_PATTERNS),
        )

    def test_book_rejects_meson_test_split_by_spaces(self):
        self.assertIn(
            "documentation test command",
            labels("meson   test -C bbdir", PROCESS_PATTERNS),
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


SPAN = (
    '<code class="docutils literal notranslate">'
    '<span class="pre">{first}</span> <span class="pre">{second}</span>'
    "</code>"
)


def _probe(root: Path, body: str) -> None:
    build = root / "docs" / "build"
    build.mkdir(parents=True)
    (build / "probe.html").write_text("<p>" + body + "</p>\n", encoding="utf-8")


def _main_stderr(root: Path) -> tuple[int, str]:
    buf = io.StringIO()
    with contextlib.redirect_stderr(buf):
        code = main(root)
    return code, buf.getvalue()


class ReaderTextEntry(unittest.TestCase):
    def test_main_rejects_words_split_across_literal_spans(self):
        """Sphinx puts one span per word. The scan has to join them.

        Dropping ``normalised_text`` from ``main`` leaves the span gaps
        as several spaces, and this phrase no longer matches.
        """
        spans = " ".join(
            f'<span class="pre">{word}</span>'
            for word in ("out", "of", "the", "repository")
        )
        html = (
            '<code class="docutils literal notranslate">' + spans + "</code>"
        )
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            _probe(root, html)
            code, err = _main_stderr(root)
        self.assertEqual(code, 1)
        self.assertIn("repository policy", err)
        self.assertIn("out of the repository", err)

    def test_main_rejects_meson_test_wrapped_across_a_line(self):
        """A literal that breaks between the two words is still the command.

        ``\\bmeson test\\b`` misses it, because the visible text keeps the
        newline between the spans.
        """
        html = SPAN.format(first="meson", second="test").replace(
            "</span> <span", "</span>\n<span"
        )
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            _probe(root, html)
            code, err = _main_stderr(root)
        self.assertEqual(code, 1)
        self.assertIn("documentation test command", err)
        self.assertIn("meson", err)
        self.assertIn("test", err)


if __name__ == "__main__":
    unittest.main()
