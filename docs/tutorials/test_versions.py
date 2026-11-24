"""The release pin is one string. These files have to agree."""

import ast
import importlib.util
import json
import re
from pathlib import Path


def _field(path: Path, pattern: str) -> str:
    text = path.read_text(encoding="utf-8")
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
    assert isinstance(value, ast.Call), ast.dump(value)
    assert not value.args and not value.keywords, ast.dump(value)
    if isinstance(value.func, ast.Name):
        module = imported[value.func.id]
        spec = importlib.util.spec_from_file_location(module, root / "scripts" / (module + ".py"))
        loaded = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(loaded)
        return getattr(loaded, value.func.id)()
    assert (
        isinstance(value.func, ast.Attribute)
        and isinstance(value.func.value, ast.Name)
    ), ast.dump(value)
    module_name = value.func.value.id
    spec = importlib.util.spec_from_file_location(module_name, root / "docs" / (module_name + ".py"))
    loaded = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(loaded)
    return getattr(loaded, value.func.attr)()


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
        "sphinx": _sphinx_release(root),
        "changelog": _field(root / "CHANGELOG.md", r"^## \[([^\]]+)\] - \d{4}-\d{2}-\d{2}$"),
        "release-notes": _field(
            root / "docs/orgmode/reference/release-2.14.org", r"^#\+TITLE:\s*(\S+) release notes$"
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
        assert not re.search(r"^[ \t]*PROJECT_NUMBER\s*=", text, re.M), rel
