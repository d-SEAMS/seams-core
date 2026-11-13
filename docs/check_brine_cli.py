#!/usr/bin/env python3
"""The paper brine CLI writes a png at the requested millimetre width."""

import sys
import tempfile
from pathlib import Path

from PIL import Image

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "docs" / "examples"))

import render_brine  # noqa: E402

CACHE = (
    ROOT / "docs" / ".cache" / "brine-front-230K"
    / "4828e33ecf419c70643b2d466bfa4f283b4138dc" / "T230_m20_r0_tend.xyz"
)


def main():
    if not CACHE.is_file():
        raise SystemExit(f"brine frame is missing: {CACHE}")
    width_mm = 40.0
    expected = render_brine.panel_width_px(width_mm)
    with tempfile.TemporaryDirectory() as tmp:
        out = Path(tmp) / "brine.png"
        render_brine.main(["--xyz", str(CACHE), "--out", str(out), "--width-mm", str(width_mm)])
        image = Image.open(out)
        if image.width < expected or image.height < expected:
            raise SystemExit(f"{image.size} is smaller than the {expected} px panel")
    print(f"brine cli wrote a figure wider than the {expected} px panel")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
