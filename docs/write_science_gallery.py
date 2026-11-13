#!/usr/bin/env python3
"""Write the science-figure page from the shared census CSV."""

import csv
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
CSV = ROOT / "paper" / "tex" / "figures" / "data" / "science" / "figure-census.csv"
PAGE = ROOT / "docs" / "source" / "science" / "gallery.rst"


def load(path=CSV):
    rows = []
    with path.open(newline="") as handle:
        for row in csv.DictReader(handle):
            rows.append((row["figure"], row["quantity"], int(row["count"])))
    return rows


def _count(rows, figure, quantity):
    for name, key, count in rows:
        if name == figure and key == quantity:
            return count
    raise KeyError((figure, quantity))


def page_text(rows):
    interface = (
        f"{_count(rows, 'interface', 'hexagonal')} hexagonal, "
        f"{_count(rows, 'interface', 'water')} liquid, "
        f"{_count(rows, 'interface', 'interfacial')} interfacial, "
        f"{_count(rows, 'interface', 'clathrate')} clathrate, "
        f"{_count(rows, 'interface', 'interClathrate')} interfacial clathrate, "
        f"and {_count(rows, 'interface', 'cubic')} cubic."
    )
    brine = (
        f"{_count(rows, 'brine', 'hexagonal')} hexagonal, "
        f"{_count(rows, 'brine', 'water')} liquid, "
        f"{_count(rows, 'brine', 'interfacial')} interfacial, "
        f"{_count(rows, 'brine', 'clathrate')} clathrate, "
        f"{_count(rows, 'brine', 'interClathrate')} interfacial clathrate, "
        f"and {_count(rows, 'brine', 'cubic')} cubic."
    )
    sii = (
        f"{_count(rows, 'sii', 'waters')} waters, "
        f"{_count(rows, 'sii', '512')} of :math:`5^{{12}}`, "
        f"and {_count(rows, 'sii', '51264')} of :math:`5^{{12}}6^{{4}}`."
    )
    sh = (
        f"{_count(rows, 'sh', 'waters')} waters, "
        f"{_count(rows, 'sh', '512')} of :math:`5^{{12}}`, "
        f"{_count(rows, 'sh', 'sh')} of :math:`4^{{3}}5^{{6}}6^{{3}}`, "
        f"and {_count(rows, 'sh', '51268')} of :math:`5^{{12}}6^{{8}}`."
    )
    return f"""Science figures
===============

The four renders used by the science tutorials. Each caption is written
from ``paper/tex/figures/data/science/figure-census.csv``. The sII cage
counts are the same numbers as ``dsv@w@sII`` in
``paper/tex/figures/data/macros.tex``.

.. figure:: figures/interface.png
   :alt: Interface slab. {interface}

   Interface. {interface}

.. figure:: figures/brine.png
   :alt: Brine front. {brine}

   Brine front. {brine}

.. figure:: figures/sh.png
   :alt: sH cages. {sh}

   sH. {sh}

.. figure:: figures/sii.png
   :alt: sII cages. {sii}

   sII. {sii}
"""


def main(argv):
    text = page_text(load())
    check = "--check" in argv
    if check:
        if not PAGE.is_file() or PAGE.read_text() != text:
            print("science gallery does not match figure-census.csv", file=sys.stderr)
            return 1
        return 0
    PAGE.write_text(text)
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
