"""Unit tests for the cross-language workflow normalizer."""

from __future__ import annotations

import importlib.util
from pathlib import Path
import unittest
from unittest import mock


SCRIPT = Path(__file__).with_name("workflow_parity.py")
SPEC = importlib.util.spec_from_file_location("workflow_parity", SCRIPT)
if SPEC is None or SPEC.loader is None:
    raise RuntimeError(f"cannot load {SCRIPT}")
workflow_parity = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(workflow_parity)


class CliParserTests(unittest.TestCase):
    def test_cage_collection_uses_short_k_option(self) -> None:
        def fake_command(_binary, command, _path, *_options):
            if command == "read":
                return {"read.nop": 1}
            return {}

        with mock.patch.object(
            workflow_parity, "_cli_command", side_effect=fake_command
        ) as command:
            workflow_parity.collect_cli(
                Path("seams"), Path("water"), Path("ice"), Path("ions")
            )

        cage_args = next(
            call.args for call in command.call_args_list if call.args[1] == "cages"
        )
        self.assertIn("-k", cage_args)
        self.assertNotIn("--k", cage_args)

    def test_scalar_commands(self) -> None:
        cases = {
            "read": ("nop 250 frame 1 box 40 40 180.203", {"read.nop": 250}),
            "chill-plus": (
                "nop 250 interClathrate 12 water 238",
                {
                    "chill_plus.nop": 250,
                    "chill_plus.interClathrate": 12,
                    "chill_plus.water": 238,
                },
            ),
            "cages": (
                "nop 4096 graph seeded hexagonal 0 cubic 4096 water 0",
                {
                    "cages.nop": 4096,
                    "cages.ih": 0,
                    "cages.ic": 4096,
                    "cages.water": 0,
                },
            ),
            "cn": (
                "# site-site\n# types 2 2 cutoff 6 cn 3.75 nI 250 nJ 250 volume 1",
                {"cn.value": 3.75},
            ),
            "hbonds": ("nop 250 hbonds 312", {"hbonds.nop": 250, "hbonds.count": 312}),
            "pairs": (
                "contact-pair count 2 nCation 2 nAnion 2",
                {"pairs.count": 2, "pairs.n_cation": 2, "pairs.n_anion": 2},
            ),
            "domains": (
                "subset polar n 2 largest 1 P 0.5",
                {"domain.n": 2, "domain.largest": 1, "domain.percolation": 0.5},
            ),
        }
        for command, (text, expected) in cases.items():
            with self.subTest(command=command):
                parsed = workflow_parity.parse_cli_text(command, text)
                for key, value in expected.items():
                    self.assertEqual(parsed[key], value)

    def test_histogram_commands(self) -> None:
        rdf = workflow_parity.parse_cli_text(
            "rdf",
            "# r g count\n# types 2 2 rmax 1 bins 2 volume 10\n0.25 0 0\n0.75 1.5 4\n",
        )
        self.assertEqual(rdf["rdf.r"], [0.25, 0.75])
        self.assertEqual(rdf["rdf.g"], [0.0, 1.5])

        density = workflow_parity.parse_cli_text(
            "density-z", "# z rho\n2.5 0.004\n7.5 0\n"
        )
        self.assertEqual(density["density.centres"], [2.5, 7.5])
        self.assertEqual(density["density.rho"], [0.004, 0.0])


class ComparisonTests(unittest.TestCase):
    def test_numeric_tolerance_and_mismatch_reporting(self) -> None:
        reference = {"count": 2, "values": [0.25, 1.5]}
        close = {"count": 2, "values": [0.25000001, 1.500001]}
        checks = workflow_parity.compare_interfaces(
            {"cli": reference, "python": close, "lua": reference}
        )
        self.assertTrue(all(check["pass"] for check in checks))

        checks = workflow_parity.compare_interfaces(
            {
                "cli": reference,
                "python": {"count": 3, "values": [0.25, 1.5]},
                "lua": reference,
            }
        )
        failed = [check for check in checks if not check["pass"]]
        self.assertEqual([check["key"] for check in failed], ["count"])


if __name__ == "__main__":
    unittest.main()
