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
    dataset_revision,
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

    def test_rejects_dataset_revision_in_rst_and_org_markup(self):
        rev = dataset_revision()
        self.assertIsNotNone(rev)
        old_samples = (
            (
                "The helper uses ``d-SEAMS/seams-core-inputs`` revision\n"
                f"``{rev}`` and verifies the pinned SHA-256 digest."
            ),
            (
                "The helper uses =d-SEAMS/seams-core-inputs= revision\n"
                f"={rev}= and verifies the pinned SHA-256 digest."
            ),
            (
                "Large tutorial and benchmark trajectories are in the public\n"
                "~d-SEAMS/seams-core-inputs~ Hugging Face dataset at revision\n"
                f"~{rev}~. Run ~python3 input/fetch.py~."
            ),
            f"Download from ``d-SEAMS/seams-core-inputs`` (``{rev}``).",
        )
        for text in old_samples:
            with self.subTest(text=text):
                found = labels(text, PROCESS_PATTERNS)
                self.assertIn("dataset revision pin", found)
                if "revision" in text:
                    self.assertIn("repository provenance", found)

        new_samples = (
            (
                "The helper downloads the trajectory from\n"
                "``d-SEAMS/seams-core-inputs`` and verifies its SHA-256 digest."
            ),
            (
                "Large tutorial and benchmark trajectories live in the public\n"
                "~d-SEAMS/seams-core-inputs~ Hugging Face dataset. Run\n"
                "~python3 input/fetch.py~ from the checkout to download all three files\n"
                "and verify their SHA-256 digests."
            ),
        )
        for text in new_samples:
            with self.subTest(text=text):
                self.assertFalse(labels(text, PROCESS_PATTERNS))

    def test_rejects_any_rst_double_backtick_revision_hash(self):
        old_rst = "The helper uses dataset revision\n``deadbeef1234567890abcdef``."
        self.assertIn("repository provenance", labels(old_rst, PROCESS_PATTERNS))
        old_org = "The helper uses dataset revision\n=deadbeef1234567890abcdef=."
        self.assertIn("repository provenance", labels(old_org, PROCESS_PATTERNS))


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

    def test_main_rejects_dataset_revision_in_rst_org_and_html(self):
        rev = dataset_revision()
        self.assertIsNotNone(rev)
        old_rst = (
            "The driver invokes ``input/fetch.py``.\n"
            "The helper uses ``d-SEAMS/seams-core-inputs`` revision\n"
            f"``{rev}`` and verifies the SHA-256 digest.\n"
        )
        new_rst = (
            "The driver invokes ``input/fetch.py``.\n"
            "The helper downloads the trajectory from ``d-SEAMS/seams-core-inputs``\n"
            "and verifies its SHA-256 digest.\n"
        )
        old_org = (
            "The helper uses =d-SEAMS/seams-core-inputs= revision\n"
            f"={rev}= and verifies the SHA-256 digest.\n"
        )
        new_org = (
            "The helper downloads the trajectory from =d-SEAMS/seams-core-inputs=\n"
            "and verifies its SHA-256 digest.\n"
        )
        old_html = (
            'The helper uses <code class="docutils literal notranslate">'
            '<span class="pre">d-SEAMS/seams-core-inputs</span></code> revision\n'
            '<code class="docutils literal notranslate">'
            f'<span class="pre">{rev}</span></code> and verifies the SHA-256 digest.'
        )
        new_html = (
            'The helper downloads the trajectory from '
            '<code class="docutils literal notranslate">'
            '<span class="pre">d-SEAMS/seams-core-inputs</span></code> '
            'and verifies its SHA-256 digest.'
        )

        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            src = root / "docs" / "source" / "explanation"
            src.mkdir(parents=True)
            (src / "rings-one-thread.rst").write_text(old_rst, encoding="utf-8")
            code, err = _main_stderr(root)
            self.assertEqual(code, 1)
            self.assertIn("repository provenance", err)
            self.assertIn("dataset revision pin", err)

            (src / "rings-one-thread.rst").write_text(new_rst, encoding="utf-8")
            code, err = _main_stderr(root)
            self.assertEqual(code, 0, err)

        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            org = root / "docs" / "orgmode" / "explanation"
            org.mkdir(parents=True)
            (org / "rings-one-thread.org").write_text(old_org, encoding="utf-8")
            code, err = _main_stderr(root)
            self.assertEqual(code, 1)
            self.assertIn("repository provenance", err)
            self.assertIn("dataset revision pin", err)

            (org / "rings-one-thread.org").write_text(new_org, encoding="utf-8")
            code, err = _main_stderr(root)
            self.assertEqual(code, 0, err)

        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            _probe(root, old_html)
            code, err = _main_stderr(root)
            self.assertEqual(code, 1)
            self.assertIn("repository provenance", err)
            self.assertIn("dataset revision pin", err)

            (root / "docs" / "build" / "probe.html").write_text(
                "<p>" + new_html + "</p>\n", encoding="utf-8"
            )
            code, err = _main_stderr(root)
            self.assertEqual(code, 0, err)

    def test_main_reads_dataset_revision_from_root_fetch_py(self):
        custom_rev = "0123456789abcdef0123456789abcdef01234567"
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            (root / "input").mkdir(parents=True)
            (root / "input" / "fetch.py").write_text(
                f'REVISION = "{custom_rev}"\n', encoding="utf-8"
            )
            _probe(
                root,
                f'Dataset pin <code class="docutils literal notranslate">'
                f'<span class="pre">{custom_rev}</span></code>.',
            )
            code, err = _main_stderr(root)
            self.assertEqual(code, 1)
            self.assertIn("dataset revision pin", err)


if __name__ == "__main__":
    unittest.main()
