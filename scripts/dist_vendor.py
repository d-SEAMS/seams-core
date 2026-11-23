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


def fetch_pin(url, pin, dest):
    """Fetch ``pin`` and leave ``dest`` checked out at that commit.

    ``git rev-parse HEAD`` must equal ``pin``. Anything else is a failed dist.
    """
    dest.mkdir(parents=True)
    subprocess.run(["git", "init", "--quiet", str(dest)], check=True)
    subprocess.run(["git", "-C", str(dest), "remote", "add", "origin", url], check=True)
    subprocess.run(
        ["git", "-C", str(dest), "fetch", "--quiet", "--depth", "1", "origin", pin],
        check=True,
    )
    subprocess.run(
        ["git", "-C", str(dest), "checkout", "--quiet", "--detach", "FETCH_HEAD"],
        check=True,
    )
    got = subprocess.run(
        ["git", "-C", str(dest), "rev-parse", "HEAD"],
        check=True,
        capture_output=True,
        text=True,
    ).stdout.strip()
    if got != pin:
        raise SystemExit("{} checked out {} but the pin is {}".format(url, got, pin))


def place_binding(dest, staging, rename=None):
    """Move a fetched binding onto ``dest`` without a missing-dest window."""
    sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "docs"))
    from binding_langs import install_snapshot

    install_snapshot(dest, staging, rename=rename)


def vendor_bindings(source_root, dist_root):
    """Put the pinned language-binding commits into the dist cache.

    docs/check_binding_api.py reads those trees. The offline image runs
    with SEAMS_NO_NETWORK=1, so the clones have to already be in the tarball.
    Each tree records the pin in ``.seams-vendored-rev`` before ``.git`` goes,
    and the checkout is the pinned commit rather than whatever HEAD is today.
    """
    sys.path.insert(0, str(source_root / "docs"))
    from binding_langs import PINS, REPOS, VENDORED_REV

    if set(PINS) != set(REPOS):
        raise SystemExit("binding pins and repositories do not name the same set")
    dest_root = dist_root / "docs" / ".cache" / "bindings"
    dest_root.mkdir(parents=True, exist_ok=True)
    for name, url in REPOS.items():
        pin = PINS[name]
        if len(pin) != 40 or any(ch not in "0123456789abcdef" for ch in pin):
            raise SystemExit("{} pin is not a 40-digit SHA".format(name))
        staging = dest_root / (".{}.partial".format(name))
        if staging.exists():
            shutil.rmtree(staging)
        try:
            fetch_pin(url, pin, staging)
            (staging / VENDORED_REV).write_text(pin + "\n", encoding="utf-8")
            git_dir = staging / ".git"
            if git_dir.is_dir():
                shutil.rmtree(git_dir)
            elif git_dir.exists():
                git_dir.unlink()
        except BaseException:
            shutil.rmtree(staging, ignore_errors=True)
            raise
        dest = dest_root / name
        # Same crash-safe replace as the online checkout. Deleting dest
        # before the rename left a killed dist with no binding tree.
        place_binding(dest, staging)


def write_tracked_manifest(source_root, dist_root):
    """Record tracked files for checks that run from an unpacked source tarball."""
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    from tracked_files import git_env

    listed = subprocess.run(
        ["git", "-C", str(source_root), "ls-files", "-z"],
        check=True,
        capture_output=True,
        env=git_env(),
    )
    out = dist_root / "docs" / ".cache" / "tracked-files.zlist"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_bytes(listed.stdout)


def main():
    source_root = Path(os.environ["MESON_SOURCE_ROOT"])
    dist_root = Path(os.environ["MESON_DIST_ROOT"])
    write_tracked_manifest(source_root, dist_root)
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
