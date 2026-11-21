#!/usr/bin/env python3
# Copyright (c) d-SEAMS developers
"""Put downloaded build inputs and test data in a Meson dist tree."""

import hashlib
import os
import shutil
import subprocess
import sys
import urllib.request
from pathlib import Path

CHEMFILES = (
    "https://github.com/chemfiles/chemfiles/archive/refs/tags/0.10.4.tar.gz",
    "b8232ddaae2953538274982838aa6c2df87d300f7e2f80e92c171581e06325ba",
    "cmake/vendor/chemfiles-0.10.4.tar.gz",
)
BRINE_REV = "4828e33ecf419c70643b2d466bfa4f283b4138dc"
DATA = (
    (
        "https://huggingface.co/datasets/d-SEAMS/brine-front-230K/resolve/"
        + BRINE_REV
        + "/data/frames/T230_m20_r0_t0.xyz",
        "dc55832fcae664aa5e528283b3be88762efdf1dabf4db1520bb1c806401b6128",
        "docs/.cache/brine-front-230K/" + BRINE_REV + "/T230_m20_r0_t0.xyz",
    ),
    (
        "https://huggingface.co/datasets/d-SEAMS/brine-front-230K/resolve/"
        + BRINE_REV
        + "/data/frames/T230_m20_r0_tend.xyz",
        "dcbc05f81f321b45195d259cd4c22e69340f9c1ae57d93c47a74c995c58dd05d",
        "docs/.cache/brine-front-230K/" + BRINE_REV + "/T230_m20_r0_tend.xyz",
    ),
    (
        "https://ndownloader.figshare.com/files/20418309",
        "c3c75050bb7762cf0d3417e2fbab326422f23ea04a64b6f69eaf33334e4d1403",
        "docs/.cache/figshare-11448702/nucleation.lammpstrj",
    ),
)
VENDOR_CONFIG = """\
[source.crates-io]
replace-with = "vendored-sources"

[source.vendored-sources]
directory = "vendor"
"""


def ignore_checkout(_directory, names):
    return [name for name in names if name in {".git", "target"}]


def copy_subprojects(source_root, dist_root):
    source = source_root / "subprojects"
    dest = dist_root / "subprojects"
    dest.mkdir(exist_ok=True)
    for child in sorted(source.iterdir()):
        if not child.is_dir() or child.name in {"packagecache", "packagefiles"}:
            continue
        target = dest / child.name
        if target.exists():
            shutil.rmtree(target)
        shutil.copytree(child, target, ignore=ignore_checkout)


def drop_build_trees(path):
    for name in ("target", ".git"):
        for found in sorted(path.rglob(name), key=lambda item: len(item.parts), reverse=True):
            if found.is_dir():
                shutil.rmtree(found)
            elif found.exists():
                found.unlink()


def vendor_cargo(project):
    if not (project / "Cargo.toml").is_file():
        return
    cargo = shutil.which("cargo")
    if cargo is None:
        raise SystemExit("cargo is required to vendor Rust subprojects")
    subprocess.run([cargo, "vendor", "vendor"], cwd=project, check=True)
    cargo_dir = project / ".cargo"
    cargo_dir.mkdir(exist_ok=True)
    config = cargo_dir / "config.toml"
    existing = config.read_text() if config.exists() else ""
    if "vendored-sources" not in existing:
        if existing and not existing.endswith("\n"):
            existing += "\n"
        config.write_text(existing + VENDOR_CONFIG)


def fetch_checked(dist_root, item):
    url, expected, relative = item
    dest = dist_root / relative
    dest.parent.mkdir(parents=True, exist_ok=True)
    got = hashlib.sha256(dest.read_bytes()).hexdigest() if dest.is_file() else ""
    if got != expected:
        part = dest.with_name(dest.name + ".part")
        urllib.request.urlretrieve(url, part)
        got = hashlib.sha256(part.read_bytes()).hexdigest()
        if got != expected:
            part.unlink(missing_ok=True)
            raise SystemExit("{} has sha256 {}".format(url, got))
        part.replace(dest)


def vendor_bindings(source_root, dist_root):
    """Shallow-clone the language bindings into the dist cache.

    docs/check_binding_api.py reads those trees. The offline image runs
    with SEAMS_NO_NETWORK=1, so the clones have to already be in the tarball.
    """
    sys.path.insert(0, str(source_root / "docs"))
    from binding_langs import REPOS

    dest_root = dist_root / "docs" / ".cache" / "bindings"
    dest_root.mkdir(parents=True, exist_ok=True)
    for name, url in REPOS.items():
        dest = dest_root / name
        if dest.exists():
            shutil.rmtree(dest)
        subprocess.run(
            ["git", "clone", "--quiet", "--depth", "1", url, str(dest)],
            check=True,
        )
        git_dir = dest / ".git"
        if git_dir.is_dir():
            shutil.rmtree(git_dir)


def main():
    source_root = Path(os.environ["MESON_SOURCE_ROOT"])
    dist_root = Path(os.environ["MESON_DIST_ROOT"])
    copy_subprojects(source_root, dist_root)
    subprojects = dist_root / "subprojects"
    drop_build_trees(subprojects)
    for child in sorted(subprojects.iterdir()):
        if child.is_dir():
            vendor_cargo(child)
    for item in (CHEMFILES,) + DATA:
        fetch_checked(dist_root, item)
    vendor_bindings(source_root, dist_root)


if __name__ == "__main__":
    main()
