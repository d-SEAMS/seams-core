#!/usr/bin/env python3
# Copyright (c) d-SEAMS developers
"""Prove that bumping meson.build's project version propagates to every check."""

from __future__ import annotations

import ast
import importlib.util
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

BUMPED_VERSION = "9.8.7"
BUMPED_MAJOR, BUMPED_MINOR, BUMPED_PATCH = 9, 8, 7
BUMPED_PACKED = (BUMPED_MAJOR << 16) | (BUMPED_MINOR << 8) | BUMPED_PATCH

COPY_PATHS = (
    "meson.build",
    "CHANGELOG.md",
    "CITATION.cff",
    "codemeta.json",
    "pixi.toml",
    "context7.json",
    "AGENTS.md",
    "README.md",
    "llms.txt",
    "Doxyfile-prj.cfg",
    "Doxyfile.cfg",
    "recipe",
    "skills",
    "packaging",
    "nix",
    "scripts",
    "include",
    "bindings",
    "docs",
)


def _copy_tree(src_root: Path, dst_root: Path) -> None:
    ignore = shutil.ignore_patterns(
        "__pycache__",
        "build",
        "doctrees",
        "doxygen-html",
        "doxygen-include",
        "lang-cpp",
        "linkcheck",
        "doctest",
        "xml",
        ".cache",
    )
    for rel in COPY_PATHS:
        src = src_root / rel
        dst = dst_root / rel
        if src.is_dir():
            shutil.copytree(src, dst, ignore=ignore)
        else:
            dst.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(src, dst)


def _load_module(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, str(path))
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _meson_project_version(root: Path) -> str:
    text = (root / "meson.build").read_text(encoding="utf-8")
    match = re.search(
        r"^project\([^)]*?version:\s*'([0-9]+\.[0-9]+\.[0-9]+)'",
        text,
        re.M,
    )
    assert match, "meson.build has no project version"
    return match.group(1)


def _bump_meson(root: Path, version: str) -> None:
    meson_path = root / "meson.build"
    text = meson_path.read_text(encoding="utf-8")
    updated, count = re.subn(
        r"(^project\([^)]*?version:\s*')[^']+(')",
        r"\g<1>" + version + r"\g<2>",
        text,
        count=1,
        flags=re.M,
    )
    assert count == 1, "failed to bump meson.build version"
    meson_path.write_text(updated, encoding="utf-8")


def _sync_release_files(root: Path, version: str) -> None:
    pixi = root / "pixi.toml"
    pixi.write_text(
        re.sub(r'^(version\s*=\s*)"[^"]+"', r'\1"{}"'.format(version), pixi.read_text(encoding="utf-8"), count=1, flags=re.M),
        encoding="utf-8",
    )
    recipe = root / "recipe" / "recipe.yaml"
    recipe.write_text(
        re.sub(r'^(\s+version:\s*)"[^"]+"', r'\1"{}"'.format(version), recipe.read_text(encoding="utf-8"), count=1, flags=re.M),
        encoding="utf-8",
    )
    cff = root / "CITATION.cff"
    cff.write_text(
        re.sub(r'^(version:\s*)"[^"]+"', r'\1"{}"'.format(version), cff.read_text(encoding="utf-8"), count=1, flags=re.M),
        encoding="utf-8",
    )
    codemeta_path = root / "codemeta.json"
    codemeta = json.loads(codemeta_path.read_text(encoding="utf-8"))
    codemeta["version"] = version
    codemeta_path.write_text(json.dumps(codemeta, indent=2) + "\n", encoding="utf-8")

    skill = root / "skills" / "d-seams" / "SKILL.md"
    skill.write_text(
        re.sub(r'^(\s+version:\s*)"[^"]+"', r'\1"{}"'.format(version), skill.read_text(encoding="utf-8"), count=1, flags=re.M),
        encoding="utf-8",
    )
    ygg = root / "packaging" / "yggdrasil" / "L" / "libseams" / "build_tarballs.jl"
    ygg.write_text(
        re.sub(r'^(version\s*=\s*v)"[^"]+"', r'\1"{}"'.format(version), ygg.read_text(encoding="utf-8"), count=1, flags=re.M),
        encoding="utf-8",
    )
    notes = root / "docs/orgmode/reference/release-2.14.org"
    notes.write_text(
        re.sub(
            r"^(#\+TITLE:\s*)\S+(\s+release notes)$",
            r"\g<1>" + version + r"\g<2>",
            notes.read_text(encoding="utf-8"),
            count=1,
            flags=re.M,
        ),
        encoding="utf-8",
    )
    changelog = root / "CHANGELOG.md"
    released = re.search(r'^date-released:\s*"([^"]+)"', cff.read_text(encoding="utf-8"), re.M).group(1)
    changelog.write_text(
        "## [{}] - {}\n\n- Bump test.\n\n".format(version, released) + changelog.read_text(encoding="utf-8"),
        encoding="utf-8",
    )


def _test_release_sync(tmp: Path) -> None:
    script = str(tmp / "scripts" / "check_release_version.py")
    drift = subprocess.run([sys.executable, script], capture_output=True, text=True)
    assert drift.returncode != 0, "check_release_version.py must fail before metadata files are synced"
    for expected in (
        "pixi.toml",
        "recipe/recipe.yaml",
        "CITATION.cff",
        "codemeta.json",
        "skills/d-seams/SKILL.md",
        "packaging/yggdrasil/L/libseams/build_tarballs.jl",
        "CHANGELOG.md",
    ):
        assert expected in drift.stderr, "missing drift error for {}: {}".format(expected, drift.stderr)

    _sync_release_files(tmp, BUMPED_VERSION)
    ok = subprocess.run([sys.executable, script], capture_output=True, text=True)
    assert ok.returncode == 0, ok.stderr
    assert "release version {} agrees in 8 files".format(BUMPED_VERSION) in ok.stdout

    printed = subprocess.run([sys.executable, script, "--print"], capture_output=True, text=True)
    assert printed.returncode == 0, printed.stderr
    assert printed.stdout.strip() == BUMPED_VERSION


def _test_doxygen_generation(tmp: Path, released: str) -> None:
    subprocess.run(
        [sys.executable, str(tmp / "scripts" / "doxygen_version_header.py")],
        check=True,
        capture_output=True,
        text=True,
    )
    inc = tmp / "docs" / "source" / "doxygen-include"
    header = (inc / "seams_version.h").read_text(encoding="utf-8")
    assert "#define SEAMS_VERSION_MAJOR {}".format(BUMPED_MAJOR) in header
    assert "#define SEAMS_VERSION_MINOR {}".format(BUMPED_MINOR) in header
    assert "#define SEAMS_VERSION_PATCH {}".format(BUMPED_PATCH) in header
    assert '#define SEAMS_VERSION_STRING "{}"'.format(BUMPED_VERSION) in header
    assert "#define SEAMS_VERSION_NUMBER {}".format(BUMPED_PACKED) in header

    v_cfg = (inc / "Doxyfile_version.cfg").read_text(encoding="utf-8")
    p_cfg = (inc / "Doxyfile_version_plain.cfg").read_text(encoding="utf-8")
    assert 'PROJECT_NUMBER          = "v{}"'.format(BUMPED_VERSION) in v_cfg
    assert 'PROJECT_NUMBER          = "{}"'.format(BUMPED_VERSION) in p_cfg

    for rel in (
        "docs/source/Doxyfile_doxyyoda.cfg",
        "docs/source/Doxyfile_seams.cfg",
        "Doxyfile-prj.cfg",
        "docs/Doxyfile.lang",
    ):
        text = (tmp / rel).read_text(encoding="utf-8")
        assert "@INCLUDE" in text and "Doxyfile_version" in text, rel
        assert released not in text, rel

    doxygen_bin = os.environ.get("DOXYGEN") or shutil.which("doxygen")
    if doxygen_bin:
        subprocess.run([doxygen_bin, "docs/Doxyfile.lang"], cwd=tmp, check=True, capture_output=True, text=True)
        index_html = (tmp / "docs" / "lang-cpp" / "html" / "index.html").read_text(encoding="utf-8")
        assert "v" + BUMPED_VERSION in index_html


def _test_discoverability_and_conf(tmp: Path, released: str) -> None:
    conf_text = (tmp / "docs" / "source" / "conf.py").read_text(encoding="utf-8")
    assert released not in conf_text
    code = """
import ast
import json
import sys
from pathlib import Path

root = Path({root!r})
sys.path.insert(0, str(root / "docs"))
import discoverability
import check_discoverability

assert discoverability.project_version() == {ver!r}
assert check_discoverability._version_problems() == []
probs = check_discoverability.source_problems()
assert probs == [], probs

ld = json.loads(discoverability.json_ld("index", "d-SEAMS", "Desc"))
graph = ld["@graph"]
assert graph[0]["version"] == {ver!r}
assert graph[1]["softwareVersion"] == {ver!r}
assert ("d-SEAMS " + {ver!r}) in discoverability.not_found_html()

conf_src = (root / "docs" / "source" / "conf.py").read_text(encoding="utf-8")
tree = ast.parse(conf_src)
ns = {{"discoverability": discoverability}}
for node in tree.body:
    if isinstance(node, ast.Assign) and len(node.targets) == 1 and isinstance(node.targets[0], ast.Name):
        if node.targets[0].id in ("version", "release"):
            exec(compile(ast.Module(body=[node], type_ignores=[]), "<conf>", "exec"), ns)
assert ns["version"] == {ver!r}
assert ns["release"] == {ver!r}
""".format(root=str(tmp), ver=BUMPED_VERSION)
    res = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True)
    assert res.returncode == 0, res.stderr


def _test_c_api_tour_and_video_tapes(tmp: Path, released: str) -> None:
    tour = _load_module("check_c_api_tour_tmp", tmp / "docs" / "check_c_api_tour.py")
    assert tour.meson_version(tmp) == (BUMPED_MAJOR, BUMPED_MINOR, BUMPED_PATCH)
    page = (tmp / "docs" / "api" / "c.md").read_text(encoding="utf-8")
    resolved = tour.resolve_page(page, tmp)
    expected_lines = [
        "capi seams_version {}".format(BUMPED_PACKED),
        "capi seams_version_major {}".format(BUMPED_MAJOR),
        "capi seams_version_minor {}".format(BUMPED_MINOR),
        "capi seams_version_patch {}".format(BUMPED_PATCH),
        "capi seams_version_string {}".format(BUMPED_VERSION),
    ]
    for row in expected_lines:
        assert row in resolved, "missing {}".format(row)
    for ph in ("<number>", "<major>", "<minor>", "<patch>", "<version>"):
        assert ph not in resolved, "unresolved placeholder {}".format(ph)
    stale_line = "capi seams_version_string {}".format(released)
    assert stale_line not in resolved
    assert released not in page
    conf_text = (tmp / "docs" / "source" / "conf.py").read_text(encoding="utf-8")
    assert "check_c_api_tour.resolve_page(" in conf_text

    all_doc = re.findall(r"(?m)^capi \S.*$", resolved)
    tour.check_page(all_doc, page, tmp, [], "\n".join(all_doc))
    stale_doc = [
        stale_line if r.startswith("capi seams_version_string ") else r
        for r in all_doc
    ]
    failed = False
    try:
        tour.check_page(stale_doc, page, tmp, [], "\n".join(stale_doc))
    except SystemExit:
        failed = True
    assert failed, "check_page must reject stale {} output when meson is bumped".format(released)

    tapes = _load_module("video_tapes_tmp", tmp / "docs" / "video_tapes.py")
    assert tapes.project_version(tmp) == BUMPED_VERSION
    lua_dir = tmp / "docs" / "videos" / "lua"
    tapes.check_static_folder(lua_dir, BUMPED_VERSION)
    recorded = (lua_dir / "recorded.txt").read_text(encoding="utf-8")
    assert released not in recorded
    expanded_rec = tapes.expand_version(recorded, BUMPED_VERSION)
    assert "{} cubic 8".format(BUMPED_VERSION) in expanded_rec
    assert released not in expanded_rec

    tv = _load_module("test_versions_tmp", tmp / "docs" / "tutorials" / "test_versions.py")
    tv.test_version_files_agree(tmp)
    tv.test_derived_files_read_meson(tmp)


def _test_nix_extraction(tmp: Path, released: str) -> None:
    nix_text = (tmp / "nix" / "package.nix").read_text(encoding="utf-8")
    assert "builtins.readFile ../meson.build" in nix_text
    assert "\\nproject(" in nix_text
    assert released not in nix_text
    # Mirror nix/package.nix project()-anchored extraction on the bumped meson.build,
    # including ensuring a stray version: outside project() is ignored.
    meson_text = (tmp / "meson.build").read_text(encoding="utf-8") + "\n# version: '0.0.1'\n"
    after_project = ("\n" + meson_text).split("\nproject(")[1]
    project_block = after_project.split(")")[0]
    pattern = re.compile(r"^[ \t]*version:[ \t]*'([0-9]+\.[0-9]+\.[0-9]+)',?[ \t]*$")
    matches = [
        m.group(1)
        for line in project_block.splitlines()
        if (m := pattern.fullmatch(line)) is not None
    ]
    assert matches == [BUMPED_VERSION], matches


def _test_finalize_seo(tmp: Path, released: str) -> None:
    path = tmp / "docs" / "finalize_seo.py"
    text = path.read_text(encoding="utf-8")
    assert released not in text, "finalize_seo.py still contains {}".format(released)
    values = []
    for node in ast.walk(ast.parse(text)):
        if not isinstance(node, ast.Dict):
            continue
        for key, value in zip(node.keys, node.values):
            if isinstance(key, ast.Constant) and key.value == "softwareVersion":
                values.append(value)
    assert len(values) == 1, values
    call = values[0]
    assert isinstance(call, ast.Call), ast.dump(call)
    func = call.func
    name = func.id if isinstance(func, ast.Name) else getattr(func, "attr", None)
    assert name == "project_version", ast.dump(call)
    sys.path.insert(0, str(tmp / "docs"))
    try:
        seo = _load_module("finalize_seo_tmp", path)
        payload = json.loads(seo._json_ld())
    finally:
        sys.path.remove(str(tmp / "docs"))
    versions = [
        node.get("softwareVersion")
        for node in payload["@graph"]
        if "softwareVersion" in node
    ]
    assert versions == [BUMPED_VERSION], versions


def main() -> int:
    src_root = Path(sys.argv[1]) if len(sys.argv) > 1 else Path(__file__).resolve().parents[1]
    released = _meson_project_version(src_root)
    with tempfile.TemporaryDirectory() as td:
        tmp = Path(td) / "seams-core"
        _copy_tree(src_root, tmp)
        _bump_meson(tmp, BUMPED_VERSION)
        _test_release_sync(tmp)
        _test_doxygen_generation(tmp, released)
        _test_discoverability_and_conf(tmp, released)
        _test_c_api_tour_and_video_tapes(tmp, released)
        _test_nix_extraction(tmp, released)
        _test_finalize_seo(tmp, released)
    print("version derive ok (bumped to {})".format(BUMPED_VERSION))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
