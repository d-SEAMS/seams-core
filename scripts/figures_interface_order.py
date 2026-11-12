# Copyright (c) d-SEAMS developers

"""Draw the CHILL+ interface and confined-ice figures.

The binary is ``interface_order`` from the test build. Reports are the
text it prints. SVG is written with the standard library.

    python3 scripts/figures_interface_order.py \\
        --bin bbdir-eessi/tests/interface_order \\
        --out docs/orgmode/explanation/figures
"""

from __future__ import annotations

import argparse
import os
import subprocess
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def run_report(binary: Path, args: list[str]) -> str:
    env = os.environ.copy()
    lib = binary.resolve().parents[1] / "src"
    env["LD_LIBRARY_PATH"] = str(lib) + os.pathsep + env.get("LD_LIBRARY_PATH", "")
    done = subprocess.run(
        [str(binary), *args],
        cwd=ROOT,
        env=env,
        check=True,
        capture_output=True,
        text=True,
    )
    return done.stdout


def parse_blocks(text: str) -> list[dict]:
    blocks = []
    current = None
    bins = []
    frames = []
    for line in text.splitlines():
        if line.startswith("# iface"):
            if current is not None:
                current["bins"] = bins
                current["frames"] = frames
                blocks.append(current)
            current = {"report": line.split()[-1]}
            bins = []
            frames = []
            continue
        if current is None or not line or line.startswith("bin mid") or line.startswith("frame t"):
            continue
        parts = line.split()
        if parts[0] == "bin" and len(parts) >= 8:
            bins.append(
                {
                    "mid": float(parts[1]),
                    "n": int(parts[2]),
                    "disagree": int(parts[3]),
                    "crystal_pub": float(parts[4]),
                    "solid_pub": float(parts[5]),
                    "crystal_ovi": float(parts[6]),
                    "solid_ovi": float(parts[7]),
                }
            )
            continue
        if parts[0] == "frame" and len(parts) >= 12:
            frames.append([float(x) if i else x for i, x in enumerate(parts[1:])])
            continue
        key = parts[0]
        current.setdefault(key, [])
        current[key].append(parts[1:])
    if current is not None:
        current["bins"] = bins
        current["frames"] = frames
        blocks.append(current)
    return blocks


def field(block: dict, key: str, default: str = "") -> str:
    rows = block.get(key, [])
    if not rows:
        return default
    return " ".join(rows[-1])


def series(block: dict, name: str) -> list[tuple[float, float]]:
    out = []
    for row in block["bins"]:
        if row["n"] <= 0:
            continue
        out.append((row["mid"], row[name]))
    return out


def polyline(points: list[tuple[float, float]], x0, x1, y0, y1, left, top, width, height) -> str:
    if len(points) < 2:
        return ""
    coords = []
    for x, y in points:
        px = left + (x - x0) / (x1 - x0) * width
        py = top + (1.0 - (y - y0) / (y1 - y0)) * height
        coords.append(f"{px:.2f},{py:.2f}")
    return " ".join(coords)


def axis_box(title: str, xlabel: str, ylabel: str, curves: list[tuple[str, str, list]],
             x0: float, x1: float, y0: float, y1: float) -> str:
    left, top, width, height = 70, 40, 640, 320
    parts = [
        f'<text x="390" y="24" text-anchor="middle" font-family="sans-serif" font-size="16">{title}</text>',
        f'<rect x="{left}" y="{top}" width="{width}" height="{height}" fill="none" stroke="#222"/>',
        f'<text x="390" y="390" text-anchor="middle" font-family="sans-serif" font-size="13">{xlabel}</text>',
        f'<text x="18" y="200" text-anchor="middle" font-family="sans-serif" font-size="13" transform="rotate(-90 18 200)">{ylabel}</text>',
    ]
    for frac, label in ((0.0, f"{y0:.1f}"), (0.5, f"{(y0+y1)/2:.1f}"), (1.0, f"{y1:.1f}")):
        y = top + (1.0 - frac) * height
        parts.append(f'<line x1="{left-4}" y1="{y:.1f}" x2="{left}" y2="{y:.1f}" stroke="#222"/>')
        parts.append(
            f'<text x="{left-8}" y="{y+4:.1f}" text-anchor="end" font-family="sans-serif" font-size="11">{label}</text>'
        )
    for frac in (0.0, 0.5, 1.0):
        x = left + frac * width
        tick = x0 + frac * (x1 - x0)
        parts.append(f'<line x1="{x:.1f}" y1="{top+height}" x2="{x:.1f}" y2="{top+height+4}" stroke="#222"/>')
        parts.append(
            f'<text x="{x:.1f}" y="{top+height+18}" text-anchor="middle" font-family="sans-serif" font-size="11">{tick:.0f}</text>'
        )
    legend_x = left + 12
    for name, colour, points in curves:
        path = polyline(points, x0, x1, y0, y1, left, top, width, height)
        parts.append(f'<polyline fill="none" stroke="{colour}" stroke-width="2" points="{path}"/>')
        parts.append(
            f'<line x1="{legend_x}" y1="56" x2="{legend_x+22}" y2="56" stroke="{colour}" stroke-width="2"/>'
        )
        parts.append(
            f'<text x="{legend_x+26}" y="60" font-family="sans-serif" font-size="11">{name}</text>'
        )
        legend_x += 150
    return "\n".join(parts)


def write_svg(path: Path, body: str, width: int = 780, height: int = 420) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        '<?xml version="1.0" encoding="UTF-8"?>\n'
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" '
        f'viewBox="0 0 {width} {height}">\n'
        '<rect width="100%" height="100%" fill="#fff"/>\n'
        f"{body}\n</svg>\n"
    )


def profile_figure(path: Path, block: dict, title: str, xlabel: str) -> None:
    curves = [
        ("2015 crystal", "#1f4e79", series(block, "crystal_pub")),
        ("2015 solid", "#5b9bd5", series(block, "solid_pub")),
        ("OVITO crystal", "#833c0c", series(block, "crystal_ovi")),
        ("OVITO solid", "#c65911", series(block, "solid_ovi")),
    ]
    xs = [x for curve in curves for x, _ in curve[2]]
    if not xs:
        return
    body = axis_box(title, xlabel, "fraction", curves, min(xs), max(xs), 0.0, 1.0)
    write_svg(path, body)


def confined_figure(path: Path, blocks: list[dict]) -> None:
    # One row per framework. Columns are the labels that can name it.
    rows = []
    for block in blocks:
        tag = field(block, "tag")
        n = int(field(block, "n", "1"))
        def count(key: str, name: str) -> int:
            for parts in block.get(key, []):
                if parts and parts[0] == name:
                    return int(parts[1])
            return 0
        prism = int(field(block, "prism_atoms", "0")) / n
        polygon = int(field(block, "polygon_atoms", "0")) / n
        bilayer = int(field(block, "bilayer_atoms", "0")) / n
        rings = max(prism, polygon, bilayer)
        crystal = (count("count_pub", "cubic") + count("count_pub", "hexagonal")) / n
        face = count("count_ovi", "interfacial") / n
        diamond = (int(field(block, "diamond_cubic", "0")) + int(field(block, "diamond_hex", "0"))) / n
        cna = int(field(block, "cna_crystal", "0")) / n
        ptm = int(field(block, "ptm_pass", "0")) / n
        rows.append((tag, rings, crystal, face, diamond, cna, ptm))
    left, top, row_h, col_w = 120, 50, 36, 80
    labels = ["rings", "CHILL+ crystal", "OVITO interfacial", "diamond", "CNA", "PTM"]
    parts = [
        '<text x="360" y="28" text-anchor="middle" font-family="sans-serif" font-size="16">'
        "Confined ice: fraction of molecules carrying each label</text>"
    ]
    for i, label in enumerate(labels):
        x = left + i * col_w + col_w / 2
        parts.append(
            f'<text x="{x}" y="58" text-anchor="middle" font-family="sans-serif" font-size="11">{label}</text>'
        )
    for r, row in enumerate(rows):
        y = top + 28 + r * row_h
        parts.append(
            f'<text x="{left-8}" y="{y+16}" text-anchor="end" font-family="sans-serif" font-size="13">{row[0]}</text>'
        )
        for c, value in enumerate(row[1:]):
            x = left + c * col_w
            shade = int(255 - 180 * value)
            parts.append(
                f'<rect x="{x+8}" y="{y}" width="{col_w-16}" height="24" fill="rgb({shade},{shade},255)" stroke="#222"/>'
            )
            parts.append(
                f'<text x="{x+col_w/2}" y="{y+16}" text-anchor="middle" font-family="sans-serif" font-size="11">{value:.2f}</text>'
            )
    write_svg(path, "\n".join(parts), width=720, height=top + 40 + row_h * max(len(rows), 1))


def growth_figure(path: Path, block: dict) -> None:
    frames = block.get("frames", [])
    if len(frames) < 2:
        return
    # parts after "frame": index, t_ps, n, disagree, face_pub, face_ovi,
    # solid_pub, solid_ovi, thick_solid_pub, thick_solid_ovi,
    # thick_crystal_pub, thick_crystal_ovi, ...
    def keep(row, col):
        return 0.05 < row[6] < 0.95 and row[col] > 0.0
    pub = [(row[1], row[8]) for row in frames if keep(row, 8)]
    ovi = [(row[1], row[9]) for row in frames if keep(row, 9)]
    cpub = [(row[1], row[10]) for row in frames if keep(row, 10)]
    covi = [(row[1], row[11]) for row in frames if keep(row, 11)]
    xs = [p[0] for p in pub]
    ys = [p[1] for p in pub + ovi + cpub + covi]
    curves = [
        ("2015 solid", "#5b9bd5", pub),
        ("OVITO solid", "#c65911", ovi),
        ("2015 crystal", "#1f4e79", cpub),
        ("OVITO crystal", "#833c0c", covi),
    ]
    body = axis_box(
        "Ice thickness along the slab",
        "time (ps)",
        "thickness (A)",
        curves,
        min(xs),
        max(xs),
        min(ys) - 1.0,
        max(ys) + 1.0,
    )
    write_svg(path, body)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--bin", type=Path, default=ROOT / "bbdir-eessi" / "tests" / "interface_order")
    parser.add_argument("--out", type=Path, default=ROOT / "docs" / "orgmode" / "explanation" / "figures")
    parser.add_argument("--growth", type=Path, default=None)
    parser.add_argument("--dt", type=float, default=100.0)
    args = parser.parse_args()
    binary = args.bin
    nucleus = parse_blocks(run_report(binary, ["study", "input/traj/nucleation_f1.lammpstrj", "1"]))[0]
    face = parse_blocks(run_report(binary, ["study", "input/traj/interface_f0.lammpstrj", "1"]))[0]
    confined = parse_blocks(run_report(binary, ["confined"]))
    profile_figure(
        args.out / "nucleus-profile.svg",
        nucleus,
        "mW nucleus, frame 15677000",
        "radius (A)",
    )
    profile_figure(
        args.out / "interface-profile.svg",
        face,
        "mW slab, interface_f0",
        "z (A)",
    )
    confined_figure(args.out / "confined-labels.svg", confined)
    if args.growth is not None:
        growth = parse_blocks(run_report(binary, ["growth", str(args.growth), str(args.dt)]))[0]
        growth_figure(args.out / "growth-thickness.svg", growth)
    print(f"wrote figures under {args.out}")


if __name__ == "__main__":
    main()
