#!/usr/bin/env python3
# Copyright (c) d-SEAMS developers
"""Regenerating docs/api/c.md keeps the version placeholders and is idempotent."""

import re
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

sys.dont_write_bytecode = True
import capi_docs  # noqa: E402
import check_c_api_tour  # noqa: E402

DOCS = Path(__file__).resolve().parent
ROOT = DOCS.parent
RENDER = DOCS / "render_c_api.py"
PAGE = DOCS / "api" / "c.md"


def quoted_rows(page):
    """capi lines of the Functions section, in the order the tour prints them."""
    section = page.split("\n## Functions\n", 1)[1].split(capi_docs.BEGIN, 1)[0]
    return re.findall(r"(?m)^capi \S.*$", section)


def tour_rows(page, version):
    """The capi lines the tour prints for a library built at version."""
    major, minor, patch = version
    values = {
        "<number>": (major << 16) | (minor << 8) | patch,
        "<major>": major,
        "<minor>": minor,
        "<patch>": patch,
        "<version>": "{}.{}.{}".format(major, minor, patch),
    }
    rows = []
    for row in quoted_rows(page):
        name, _, value = row[len("capi "):].partition(" ")
        if capi_docs.VERSION_PLACEHOLDERS.get(name) == value:
            row = "capi {} {}".format(name, values[value])
        rows.append(row)
    return rows


def regenerate(rows, out):
    done = subprocess.run(
        [sys.executable, str(RENDER), str(out)],
        input="\n".join(rows) + "\n",
        text=True,
        capture_output=True,
    )
    if done.returncode != 0:
        raise AssertionError(done.stderr)
    return out.read_text()


class RenderCApiSelfTest(unittest.TestCase):
    def setUp(self):
        directory = tempfile.TemporaryDirectory()
        self.addCleanup(directory.cleanup)
        self.out = Path(directory.name) / "c.md"
        self.tracked = PAGE.read_text()
        self.version = check_c_api_tour.meson_version(ROOT)
        self.rows = tour_rows(self.tracked, self.version)

    def test_tracked_page_carries_placeholders(self):
        for name, mark in capi_docs.VERSION_PLACEHOLDERS.items():
            self.assertIn("capi {} {}\n".format(name, mark), self.tracked)

    def test_regenerated_page_keeps_placeholders(self):
        text = regenerate(self.rows, self.out)
        for name, mark in capi_docs.VERSION_PLACEHOLDERS.items():
            self.assertIn("capi {} {}\n".format(name, mark), text)
        literals = [row for row in self.rows if row.split()[1] in capi_docs.VERSION_PLACEHOLDERS]
        self.assertEqual(len(literals), len(capi_docs.VERSION_PLACEHOLDERS))
        for row in literals:
            self.assertNotIn(row, text)
        self.assertNotIn("{}.{}.{}".format(*self.version), text)

    def test_regeneration_is_idempotent(self):
        first = regenerate(self.rows, self.out)
        second = regenerate(tour_rows(first, self.version), self.out)
        self.assertEqual(first, second)

    def test_regeneration_does_not_depend_on_the_version(self):
        first = regenerate(self.rows, self.out)
        bumped = regenerate(tour_rows(self.tracked, (9, 87, 65)), self.out)
        self.assertEqual(first, bumped)

    def test_regenerated_page_passes_the_tour_check(self):
        text = regenerate(self.rows, self.out)
        check_c_api_tour.check_page(self.rows, text, ROOT, [], "\n".join(self.rows))

    def test_regeneration_quotes_the_same_tour_lines(self):
        text = regenerate(self.rows, self.out)
        self.assertEqual(quoted_rows(text), quoted_rows(self.tracked))

    def test_regeneration_leaves_the_tracked_page_alone(self):
        regenerate(self.rows, self.out)
        self.assertEqual(PAGE.read_text(), self.tracked)


if __name__ == "__main__":
    unittest.main()
