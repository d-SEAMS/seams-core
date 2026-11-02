"""The release pin is one string. These files have to agree."""

import re
from pathlib import Path


def _field(path: Path, pattern: str) -> str:
    text = path.read_text()
    match = re.search(pattern, text, re.M)
    assert match, path
    return match.group(1).lstrip("v")


def test_version_files_agree(root: Path):
    pins = {
        "meson": _field(root / "meson.build", r"^\s*version:\s*'([^']+)'"),
        "pixi": _field(root / "pixi.toml", r'^version = "([^"]+)"'),
        "citation": _field(root / "CITATION.cff", r'^version: "([^"]+)"'),
        "nix": _field(root / "nix/package.nix", r'^\s*version = "([^"]+)";'),
        "doxygen": _field(
            root / "docs/source/Doxyfile_seams.cfg", r'^PROJECT_NUMBER\s+=\s+"v?([^"]+)"'
        ),
    }
    assert len(set(pins.values())) == 1, pins
