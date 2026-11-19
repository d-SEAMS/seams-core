"""Release version for cog snippets, read from meson.build.

The project version in meson.build is the one place a person writes the
version. Every copy listed in ci/version/files.txt is regenerated from it:

    pixi run --frozen -e version release-cog
    pixi run --frozen -e version version-drift

Call cogapp through ``python -m cogapp``, because cocogitto also installs
a ``cog`` executable.
"""

import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]

_PROJECT = re.compile(r"^\s*version:\s*'([^']+)'", re.MULTILINE)
_match = _PROJECT.search((ROOT / "meson.build").read_text(encoding="utf-8"))
if _match is None:
    raise RuntimeError("meson.build has no project version")

VERSION = _match.group(1)
