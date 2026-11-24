#!/usr/bin/env python3
"""Offline self-test for the live registry claims checker."""

import tempfile
import unittest
from pathlib import Path

from check_registry_claims import check_registry_pages
from registry_claims import RegistryError, claims


class RegistryClaimsSelfTest(unittest.TestCase):
    def _make_page(self, text):
        directory = tempfile.TemporaryDirectory()
        self.addCleanup(directory.cleanup)
        path = Path(directory.name) / "page.rst"
        path.write_text(text, encoding="utf-8")
        return path

    def test_live_lookup_reports_unpublished_version(self):
        page = self._make_page("Install with `pip install seams-core==9.9.9`.")

        def fake_live_versions(manager, package):
            self.assertEqual((manager, package), ("pip", "seams-core"))
            return {"2.14.1"}

        errors, files_scanned, claim_count, package_count = check_registry_pages(
            [(page, "page.rst")], resolve_versions=fake_live_versions
        )
        self.assertEqual(files_scanned, 1)
        self.assertEqual(claim_count, 1)
        self.assertEqual(package_count, 1)
        self.assertEqual(len(errors), 1)
        self.assertIn("unpublished install claim", errors[0])

    def test_live_lookup_reports_registry_failure(self):
        page = self._make_page("Install with `pip install seams-core`.")

        def fake_live_versions(manager, package):
            raise RegistryError(f"{manager} {package} unreachable")

        errors, files_scanned, claim_count, package_count = check_registry_pages(
            [(page, "page.rst")], resolve_versions=fake_live_versions
        )
        self.assertEqual(files_scanned, 1)
        self.assertEqual(claim_count, 1)
        self.assertEqual(package_count, 1)
        self.assertEqual(len(errors), 1)
        self.assertIn("live registry check failed", errors[0])

    def test_local_pip_install_is_not_a_registry_claim(self):
        self.assertEqual(list(claims("Build the wheel with `pip install .`.")), [])


if __name__ == "__main__":
    unittest.main()
