#!/usr/bin/env python3
"""The paper brine CLI writes a png at the requested millimetre width."""

import io
import sys
import tempfile
from contextlib import redirect_stdout
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
    with tempfile.TemporaryDirectory() as tmp:
        for width_mm in (85.0, 170.0):
            out = Path(tmp) / f"brine-{int(width_mm)}.png"
            buf = io.StringIO()
            with redirect_stdout(buf):
                render_brine.main(["--xyz", str(CACHE), "--out", str(out), "--width-mm", str(width_mm)])
            image = Image.open(out)
            if image.width < render_brine.panel_width_px(width_mm) or image.height < 32:
                raise SystemExit(f"{width_mm:g} mm figure is only {image.size}")
            line = buf.getvalue().strip().splitlines()[-1]
            points = float(line.split("legend", 1)[1].split()[0])
            if points < 7.0:
                raise SystemExit(f"{width_mm:g} mm legend is {points:.2f} pt")
            print(f"{width_mm:g} mm {image.size} legend {points:.2f} pt")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
