#!/usr/bin/env python3
"""Fetch Niu, Yang and Parrinello's ice-nucleation dumps.

The archive is Materials Cloud record 2020.0005/v1
(doi:10.24435/materialscloud:2020.0005/v1), file Nucleation.zip.
The published MD5 of that zip is pinned below. After extraction the
script checks the SHA-256 of Nucleation_trajectory-3, the oxygen dump
the head-to-head and the walk scripts read.

Destination defaults to <repo>/data/public-ice. Override with --dest
or SEAMS_NIU_DIR (the directory that should contain
nucleation/Ice_nucleation_trajectory/).
"""

from __future__ import annotations

import argparse
import hashlib
import os
import sys
import urllib.request
import zipfile
from pathlib import Path

URL = (
    "https://archive.materialscloud.org/api/records/kc275-mc539/files/"
    "Nucleation.zip/content"
)
DOI = "10.24435/materialscloud:2020.0005/v1"
ZIP_MD5 = "bddcd4d35b1fed3353dcc320343907ec"
MEMBER = "Ice_nucleation_trajectory/Nucleation_trajectory-3"
# SHA-256 of that member, recorded by the 31 October walk provenance.
MEMBER_SHA256 = (
    "69a91774e388e8913c5c75b634afae1de306466351bc97886e0482ccebebea56"
)


def repo_root() -> Path:
    return Path(__file__).resolve().parents[1]


def md5_file(path: Path) -> str:
    digest = hashlib.md5()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def download(url: str, dest: Path) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    tmp = dest.with_suffix(dest.suffix + ".part")
    with urllib.request.urlopen(url) as response, tmp.open("wb") as out:
        while True:
            chunk = response.read(1 << 20)
            if not chunk:
                break
            out.write(chunk)
    tmp.replace(dest)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--dest",
        type=Path,
        default=Path(
            os.environ.get("SEAMS_NIU_DIR", str(repo_root() / "data" / "public-ice"))
        ),
        help="directory that will contain nucleation/Ice_nucleation_trajectory/",
    )
    parser.add_argument("--url", default=os.environ.get("SEAMS_NIU_URL", URL))
    args = parser.parse_args()
    dest = args.dest
    if not dest.is_absolute():
        dest = (Path.cwd() / dest).resolve()
    member_out = dest / "nucleation" / MEMBER
    if member_out.is_file() and sha256_file(member_out) == MEMBER_SHA256:
        print(f"present {member_out} sha256={MEMBER_SHA256}")
        return 0

    zip_path = dest / "Nucleation.zip"
    if not zip_path.is_file() or md5_file(zip_path) != ZIP_MD5:
        print(f"fetch {args.url}", flush=True)
        download(args.url, zip_path)
    got = md5_file(zip_path)
    if got != ZIP_MD5:
        print(f"Nucleation.zip md5 {got} != pinned {ZIP_MD5}", file=sys.stderr)
        return 1
    with zipfile.ZipFile(zip_path) as zf:
        names = zf.namelist()
        match = [n for n in names if n.endswith(MEMBER) or n.endswith("/" + MEMBER.split("/")[-1])]
        exact = [n for n in names if n.replace("\\", "/").endswith(MEMBER)]
        chosen = exact[0] if exact else (match[0] if len(match) == 1 else "")
        if not chosen:
            print("Nucleation.zip has no Nucleation_trajectory-3", file=sys.stderr)
            print("\n".join(names[:40]), file=sys.stderr)
            return 1
        member_out.parent.mkdir(parents=True, exist_ok=True)
        with zf.open(chosen) as src, member_out.open("wb") as dst:
            dst.write(src.read())
    got_sha = sha256_file(member_out)
    if got_sha != MEMBER_SHA256:
        print(
            f"{member_out.name} sha256 {got_sha} != pinned {MEMBER_SHA256}",
            file=sys.stderr,
        )
        return 1
    print(
        f"ok doi={DOI} zip_md5={ZIP_MD5} member={MEMBER} sha256={got_sha} path={member_out}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
