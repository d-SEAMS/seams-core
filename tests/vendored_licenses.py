#!/usr/bin/env python3
# Copyright (c) d-SEAMS developers
"""Require each vendored third_party tree to ship its real licence.

A NOTICE.md does not stand in for the licence. argum's LICENSE must be
upstream v2.9's BSD 3-Clause text, Copyright (c) 2022, Eugene Gershnik,
the file argum.h names. The conda recipe, the Yggdrasil script and the
Nix package must record that licence. csm's GPL-2.0 text stays required.
"""

from __future__ import annotations

import hashlib
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
LICENSE_NAMES = {
    "LICENSE",
    "LICENSE.md",
    "LICENSE.MD",
    "LICENSE.txt",
    "COPYING",
    "COPYING.txt",
}
SOURCE_SUFFIXES = {".h", ".hpp", ".c", ".cpp", ".cc", ".f", ".f90", ".inc"}
# sha256 of https://github.com/gershnik/argum blob LICENSE at tag v2.9
# (commit 13615b4748272ebe4955dada30f2d4438810757b) and at master.
ARGUM_SHA256 = "b6fb259a38c4b2bc717c5772f60c7cacf06eb591c4756561d73525351cbd1abd"
ARGUM_COPYRIGHT = "Copyright (c) 2022, Eugene Gershnik"
CSM_MARK = "GNU GENERAL PUBLIC LICENSE"
RECIPE = Path("recipe/recipe.yaml")
YGG = Path("packaging/yggdrasil/L/libseams/build_tarballs.jl")
NIX = Path("nix/package.nix")


def _sources(directory: Path) -> list[Path]:
    return [
        path
        for path in directory.rglob("*")
        if path.is_file() and path.suffix in SOURCE_SUFFIXES
    ]


def problems(root: Path) -> list[str]:
    found = []
    third = root / "third_party"
    if not third.is_dir():
        return ["third_party is missing"]
    for child in sorted(path for path in third.iterdir() if path.is_dir()):
        if not _sources(child):
            continue
        licences = [
            path
            for path in child.rglob("*")
            if path.is_file() and path.name in LICENSE_NAMES
        ]
        if not licences:
            found.append(
                "{} has source and no licence file (NOTICE.md does not count)".format(
                    child.relative_to(root)
                )
            )
    argum_license = third / "argum" / "LICENSE"
    if not argum_license.is_file():
        found.append("third_party/argum/LICENSE is missing")
    else:
        digest = hashlib.sha256(argum_license.read_bytes()).hexdigest()
        if digest != ARGUM_SHA256:
            found.append(
                "third_party/argum/LICENSE sha256 {} is not upstream v2.9 {}".format(
                    digest, ARGUM_SHA256
                )
            )
        text = argum_license.read_text(encoding="utf-8", errors="replace")
        if ARGUM_COPYRIGHT not in text:
            found.append("third_party/argum/LICENSE lacks {}".format(ARGUM_COPYRIGHT))
        if "BSD 3-Clause License" not in text:
            found.append("third_party/argum/LICENSE is not the BSD 3-Clause License")
    header = third / "argum" / "argum.h"
    if not header.is_file():
        found.append("third_party/argum/argum.h is missing")
    else:
        head = "\n".join(header.read_text(encoding="utf-8", errors="replace").splitlines()[:12])
        if "Eugene Gershnik" not in head or "LICENSE file" not in head:
            found.append("argum.h no longer names Eugene Gershnik and the LICENSE file")
    csm = third / "csm" / "LICENSE.MD"
    if not csm.is_file() or CSM_MARK not in csm.read_text(encoding="utf-8", errors="replace"):
        found.append("third_party/csm/LICENSE.MD is missing its GPL-2.0 text")
    recipe = root / RECIPE
    if not recipe.is_file():
        found.append("{} is missing".format(RECIPE))
    else:
        recipe_text = recipe.read_text(encoding="utf-8")
        if "third_party/argum/LICENSE" not in recipe_text:
            found.append("{} license_file does not list third_party/argum/LICENSE".format(RECIPE))
        if "BSD-3-Clause" not in recipe_text:
            found.append("{} license does not include BSD-3-Clause".format(RECIPE))
    ygg = root / YGG
    if not ygg.is_file():
        found.append("{} is missing".format(YGG))
    else:
        ygg_text = ygg.read_text(encoding="utf-8")
        if "third_party/argum/LICENSE" not in ygg_text or "install_license" not in ygg_text:
            found.append("{} does not install third_party/argum/LICENSE".format(YGG))
    nix = root / NIX
    if not nix.is_file():
        found.append("{} is missing".format(NIX))
    else:
        if "bsd3" not in nix.read_text(encoding="utf-8"):
            found.append("{} does not record the BSD-3-Clause license".format(NIX))
    return found


def _seed(root: Path, license_bytes: bytes) -> None:
    argum = root / "third_party" / "argum"
    argum.mkdir(parents=True)
    (argum / "argum.h").write_text(
        "// Copyright 2022 Eugene Gershnik\n"
        "// license that can be found in the LICENSE file\n",
        encoding="utf-8",
    )
    (argum / "NOTICE.md").write_text("see the header\n", encoding="utf-8")
    (argum / "LICENSE").write_bytes(license_bytes)
    csm = root / "third_party" / "csm"
    csm.mkdir()
    (csm / "csm.c").write_text("int csm(void) { return 0; }\n", encoding="utf-8")
    (csm / "LICENSE.MD").write_text(CSM_MARK + "\n", encoding="utf-8")
    recipe = root / RECIPE
    recipe.parent.mkdir(parents=True)
    recipe.write_text(
        "license: MIT AND BSD-3-Clause\nlicense_file:\n  - third_party/argum/LICENSE\n",
        encoding="utf-8",
    )
    ygg = root / YGG
    ygg.parent.mkdir(parents=True)
    ygg.write_text(
        'install_license LICENSE\ncp third_party/argum/LICENSE "${WORKSPACE}/srcdir/LICENSE.argum"\n',
        encoding="utf-8",
    )
    nix = root / NIX
    nix.parent.mkdir(parents=True)
    nix.write_text("license = with lib.licenses; [ mit bsd3 ];\n", encoding="utf-8")


def self_test() -> int:
    real = problems(ROOT)
    if real:
        print("\n".join(real), file=sys.stderr)
        return 1
    license_bytes = (ROOT / "third_party" / "argum" / "LICENSE").read_bytes()
    with tempfile.TemporaryDirectory() as tmp:
        root = Path(tmp)
        _seed(root, license_bytes)
        if problems(root):
            print("a complete tree was rejected", file=sys.stderr)
            return 1
        argum = root / "third_party" / "argum"
        (argum / "LICENSE").unlink()
        missing = problems(root)
        if not any("NOTICE.md does not count" in item for item in missing):
            print("NOTICE.md was accepted in place of LICENSE", file=sys.stderr)
            return 1
        (argum / "LICENSE").write_text(
            "BSD 3-Clause License\n" + ARGUM_COPYRIGHT + "\nshortened\n",
            encoding="utf-8",
        )
        if not any("not upstream v2.9" in item for item in problems(root)):
            print("a shortened licence was accepted", file=sys.stderr)
            return 1
        (argum / "LICENSE").write_bytes(license_bytes)
        recipe = root / RECIPE
        recipe.write_text(recipe.read_text().replace("third_party/argum/LICENSE", "LICENSE"))
        if not any("license_file" in item for item in problems(root)):
            print("a recipe that omits the Argum licence was accepted", file=sys.stderr)
            return 1
        csm = root / "third_party" / "csm" / "LICENSE.MD"
        csm.write_text("placeholder\n", encoding="utf-8")
        if not any("GPL-2.0" in item for item in problems(root)):
            print("a stub csm licence was accepted", file=sys.stderr)
            return 1
    print("vendored licence self-test passed")
    return 0


def main(argv: list[str]) -> int:
    if argv == ["--self-test"]:
        return self_test()
    if argv:
        print("usage: tests/vendored_licenses.py [--self-test]", file=sys.stderr)
        return 2
    found = problems(ROOT)
    if found:
        print("\n".join(found), file=sys.stderr)
        return 1
    print("vendored licences match, including argum {}".format(ARGUM_SHA256))
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
