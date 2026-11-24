#!/usr/bin/env python3
# Copyright (c) d-SEAMS developers
"""Keep the wrap overlay's backward.hpp equal to the vendored header.

The file this repository kept at src/include/external/backward.hpp is
upstream commit 1d43bfea (backward-cpp v1.6 plus the 27 commits after
that tag). The wrap downloads the v1.6 tarball. subprojects/packagefiles
/backward-cpp/backward.hpp is what patch_directory copies over that
tarball. The offline recipes copy the same directory onto a pre-seeded
subproject, and Meson does not re-apply wrap patches once meson.build
is already there, so the header has to live in the overlay.
"""

from __future__ import annotations

import hashlib
import re
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
HEADER = Path("subprojects/packagefiles/backward-cpp/backward.hpp")
# sha256 of src/include/external/backward.hpp at 2c5d0fce, which is
# byte-identical to https://github.com/bombela/backward-cpp commit
# 1d43bfeaeaa561a4f762d3cf4add853c7b324a09.
HEADER_SHA256 = "cac5328fb08732fc9c8a3c173b126a6c913f684e08bf6e4cafb61331e4b32b8e"
# Present in 1d43bfea, absent from the v1.6 tag's backward.hpp.
POST_V16_MARKER = "#include <asm/ptrace.h>"
WRAP = Path("subprojects/backward-cpp.wrap")
COPY_SITES = (
    Path("recipe/recipe.yaml"),
    Path("nix/package.nix"),
    Path("packaging/yggdrasil/L/libseams/build_tarballs.jl"),
)


def problems(root: Path) -> list[str]:
    found = []
    header = root / HEADER
    if not header.is_file():
        found.append("{} is missing".format(HEADER))
    else:
        digest = hashlib.sha256(header.read_bytes()).hexdigest()
        if digest != HEADER_SHA256:
            found.append(
                "{} sha256 {} is not the vendored header {}".format(
                    HEADER, digest, HEADER_SHA256
                )
            )
        if POST_V16_MARKER not in header.read_text(encoding="utf-8", errors="replace"):
            found.append("{} lacks the post-v1.6 {!r} change".format(HEADER, POST_V16_MARKER))
    wrap = root / WRAP
    if not wrap.is_file():
        found.append("{} is missing".format(WRAP))
    else:
        text = wrap.read_text(encoding="utf-8")
        if "patch_directory = backward-cpp" not in text:
            found.append("{} does not apply the backward-cpp overlay".format(WRAP))
        if re.search(r"(?m)^\s*diff_files\s*=", text):
            found.append(
                "{} sets diff_files, which the offline pre-seed does not re-apply".format(WRAP)
            )
    for site in COPY_SITES:
        path = root / site
        if not path.is_file():
            found.append("{} is missing".format(site))
            continue
        text = path.read_text(encoding="utf-8")
        if "packagefiles/backward-cpp" not in text:
            found.append("{} does not copy the backward-cpp overlay".format(site))
    return found


def self_test() -> int:
    real = problems(ROOT)
    if real:
        print("\n".join(real), file=sys.stderr)
        return 1
    good = (ROOT / HEADER).read_bytes()
    with tempfile.TemporaryDirectory() as tmp:
        root = Path(tmp)
        header = root / HEADER
        header.parent.mkdir(parents=True)
        header.write_bytes(b"not the vendored header\n")
        (root / WRAP).parent.mkdir(parents=True, exist_ok=True)
        (root / WRAP).write_text("patch_directory = backward-cpp\n", encoding="utf-8")
        for site in COPY_SITES:
            path = root / site
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text("packagefiles/backward-cpp\n", encoding="utf-8")
        bad = problems(root)
        if not any(HEADER_SHA256 in item for item in bad):
            print("a replaced header was accepted", file=sys.stderr)
            return 1
        header.write_bytes(good)
        if problems(root):
            print("the vendored header was rejected", file=sys.stderr)
            return 1
        wrap = root / WRAP
        wrap.write_text(
            "patch_directory = backward-cpp\ndiff_files = backward-cpp\n",
            encoding="utf-8",
        )
        if not any("diff_files" in item for item in problems(root)):
            print("diff_files was accepted", file=sys.stderr)
            return 1
    print("backward overlay self-test passed")
    return 0


def main(argv: list[str]) -> int:
    if argv == ["--self-test"]:
        return self_test()
    if argv:
        print("usage: tests/backward_overlay.py [--self-test]", file=sys.stderr)
        return 2
    found = problems(ROOT)
    if found:
        print("\n".join(found), file=sys.stderr)
        return 1
    print("backward overlay header matches {}".format(HEADER_SHA256))
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
