"""The release pin is one string. These files have to agree."""

import ast
import importlib.util
import re
from pathlib import Path


def _field(path: Path, pattern: str) -> str:
    text = path.read_text()
    match = re.search(pattern, text, re.M)
    assert match, path
    return match.group(1).lstrip("v")


def _sphinx_release(root: Path) -> str:
    """Resolve ``release`` in conf.py through its assignments, as Sphinx does."""
    tree = ast.parse((root / "docs/source/conf.py").read_text())
    assigned = {
        node.targets[0].id: node.value
        for node in tree.body
        if isinstance(node, ast.Assign)
        and len(node.targets) == 1
        and isinstance(node.targets[0], ast.Name)
    }
    imported = {
        alias.asname or alias.name: node.module
        for node in tree.body
        if isinstance(node, ast.ImportFrom)
        for alias in node.names
    }
    value = assigned["release"]
    seen = set()
    while isinstance(value, ast.Name):
        assert value.id not in seen, value.id
        seen.add(value.id)
        value = assigned[value.id]
    if isinstance(value, ast.Constant) and isinstance(value.value, str):
        return value.value
    assert isinstance(value, ast.Call) and isinstance(value.func, ast.Name), ast.dump(value)
    assert not value.args and not value.keywords, ast.dump(value)
    module = imported[value.func.id]
    spec = importlib.util.spec_from_file_location(module, root / "scripts" / (module + ".py"))
    loaded = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(loaded)
    return getattr(loaded, value.func.id)()


def test_version_files_agree(root: Path):
    pins = {
        "meson": _field(root / "meson.build", r"^\s*version:\s*'([^']+)'"),
        "pixi": _field(root / "pixi.toml", r'^version = "([^"]+)"'),
        "citation": _field(root / "CITATION.cff", r'^version: "([^"]+)"'),
        "nix": _field(root / "nix/package.nix", r'^\s*version = "([^"]+)";'),
        "doxygen": _field(
            root / "docs/source/Doxyfile_seams.cfg", r'^PROJECT_NUMBER\s+=\s+"v?([^"]+)"'
        ),
        "doxyyoda": _field(
            root / "docs/source/Doxyfile_doxyyoda.cfg", r'^PROJECT_NUMBER\s+=\s+"v?([^"]+)"'
        ),
        "doxygen-prj": _field(root / "Doxyfile-prj.cfg", r'^PROJECT_NUMBER\s+=\s+"v?([^"]+)"'),
        "sphinx": _sphinx_release(root),
        "changelog": _field(root / "CHANGELOG.md", r"^## \[([^\]]+)\] - \d{4}-\d{2}-\d{2}$"),
        "release-notes": _field(
            root / "docs/orgmode/reference/release-2.14.org", r"^#\+TITLE:\s*(\S+) release notes$"
        ),
    }
    assert len(set(pins.values())) == 1, pins
