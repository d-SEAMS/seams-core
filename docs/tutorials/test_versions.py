"""The release pin is one string. These files have to agree."""

import json
import re
from pathlib import Path


def _field(path: Path, pattern: str) -> str:
    text = path.read_text(encoding="utf-8")
    match = re.search(pattern, text, re.M)
    assert match, path
    return match.group(1).lstrip("v")


def test_version_files_agree(root: Path):
    codemeta = json.loads((root / "codemeta.json").read_text(encoding="utf-8"))
    pins = {
        "meson": _field(root / "meson.build", r"^project\([^)]*?version:\s*'([^']+)'"),
        "pixi": _field(root / "pixi.toml", r'^version = "([^"]+)"'),
        "citation": _field(root / "CITATION.cff", r'^version: "([^"]+)"'),
        "codemeta": codemeta.get("version", ""),
        "recipe": _field(root / "recipe/recipe.yaml", r'^context:\s*\n\s+version:\s*"([^"]+)"'),
        "skill": _field(
            root / "skills/d-seams/SKILL.md",
            r'^metadata:\s*\n(?:[ \t]+[^\n]*\n)*?[ \t]+version:\s*"?([^"\s]+)"?',
        ),
        "yggdrasil": _field(
            root / "packaging/yggdrasil/L/libseams/build_tarballs.jl",
            r'^version\s*=\s*v"([^"]+)"',
        ),
    }
    assert len(set(pins.values())) == 1, pins


def test_derived_files_read_meson(root: Path):
    nix_text = (root / "nix/package.nix").read_text(encoding="utf-8")
    assert "builtins.readFile ../meson.build" in nix_text
    assert "\\nproject(" in nix_text
    assert not re.search(r'^\s*version\s*=\s*"[^"]+";', nix_text, re.M)

    conf_text = (root / "docs/source/conf.py").read_text(encoding="utf-8")
    assert "version = discoverability.project_version()" in conf_text

    for rel in (
        "docs/source/Doxyfile_doxyyoda.cfg",
        "docs/source/Doxyfile_seams.cfg",
        "Doxyfile-prj.cfg",
        "docs/Doxyfile.lang",
    ):
        text = (root / rel).read_text(encoding="utf-8")
        assert "@INCLUDE" in text and "Doxyfile_version" in text, rel
        assert not re.search(r"^PROJECT_NUMBER\s*=", text, re.M), rel
