#!/usr/bin/env python3
# Copyright (c) d-SEAMS developers
# SPDX-License-Identifier: MIT
"""Copy the C API sources and Eigen into seams-sys/vendor/.

cargo package only takes files under the crate directory, and the C++
sources live at the repository root. Run this before cargo package or
cargo publish:

    python3 bindings/rust/seams-sys/vendor.py

Eigen comes from the 3.4.0 release tarball, checked against its sha256.
Pass --eigen-tarball to use a tarball that is already on disk.
"""

import argparse
import hashlib
import pathlib
import shutil
import sys
import tarfile
import tempfile
import urllib.request

HERE = pathlib.Path(__file__).resolve().parent
ROOT = HERE.parents[2]
EIGEN_VERSION = "3.4.0"
EIGEN_URL = f"https://gitlab.com/libeigen/eigen/-/archive/{EIGEN_VERSION}/eigen-{EIGEN_VERSION}.tar.gz"
EIGEN_SHA256 = "8586084f71f9bde545ee7fa6d00288b264a2b7ac3607b974e54d13e7162c1c72"
# MPL2 for the library, BSD for AlignedBox.h, Apache for BFloat16.h.
EIGEN_NOTICES = ("COPYING.MPL2", "COPYING.BSD", "COPYING.APACHE", "COPYING.README")


def sources():
    listed = []
    for line in (HERE / "sources.txt").read_text().splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        listed.append(line.split(":", 1)[1].strip() if ":" in line else line)
    return listed


def fetch(path):
    with urllib.request.urlopen(EIGEN_URL) as response, open(path, "wb") as out:
        shutil.copyfileobj(response, out)


def unpack_eigen(tarball, dest):
    digest = hashlib.sha256(pathlib.Path(tarball).read_bytes()).hexdigest()
    if digest != EIGEN_SHA256:
        raise SystemExit(f"{tarball}: sha256 {digest}, expected {EIGEN_SHA256}")
    prefix = f"eigen-{EIGEN_VERSION}/"
    keep = [prefix + "Eigen/"] + [prefix + name for name in EIGEN_NOTICES]
    with tarfile.open(tarball) as archive:
        for member in archive.getmembers():
            if not member.isfile() or not any(member.name.startswith(k) for k in keep):
                continue
            target = dest / member.name[len(prefix):]
            target.parent.mkdir(parents=True, exist_ok=True)
            with archive.extractfile(member) as src, open(target, "wb") as out:
                shutil.copyfileobj(src, out)


def main(argv):
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--eigen-tarball", help=f"eigen-{EIGEN_VERSION}.tar.gz on disk")
    args = parser.parse_args(argv)

    vendor = HERE / "vendor"
    if vendor.exists():
        shutil.rmtree(vendor)
    copies = sources() + ["include/seams.h", "include/seams_version.h.in", "LICENSE"]
    copies += [p.relative_to(ROOT).as_posix() for p in sorted((ROOT / "src/include/internal").glob("*.hpp"))]
    copies += [p.relative_to(ROOT).as_posix() for p in sorted((ROOT / "src/include/internal").glob("*.h"))]
    for rel in copies:
        target = vendor / rel
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(ROOT / rel, target)

    with tempfile.TemporaryDirectory() as tmp:
        tarball = args.eigen_tarball
        if tarball is None:
            tarball = pathlib.Path(tmp) / f"eigen-{EIGEN_VERSION}.tar.gz"
            fetch(tarball)
        unpack_eigen(tarball, vendor / "eigen")
    count = sum(1 for p in vendor.rglob("*") if p.is_file())
    print(f"vendored {count} files into {vendor.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
