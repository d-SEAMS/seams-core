"""Release version for cog snippets, read from meson.build.

Meson's rewriter updates the project version in meson.build. Every copy
listed in ci/version/files.txt is regenerated from it:

    pixi run --frozen -e version release-cog
    pixi run --frozen -e version version-drift

Call cogapp through ``python -m cogapp``, because cocogitto also installs
a ``cog`` executable.
"""

import os
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]

_PROJECT = re.compile(r"^\s*version\s*:\s*'([^']+)'", re.MULTILINE)
_match = _PROJECT.search((ROOT / "meson.build").read_text(encoding="utf-8"))
if _match is None:
    raise RuntimeError("meson.build has no project version")

VERSION = _match.group(1)

_DATE = re.compile(r'^date-released:\s*"?(\d{4}-\d{2}-\d{2})"?', re.MULTILINE)
_date_match = _DATE.search((ROOT / "CITATION.cff").read_text(encoding="utf-8"))
if _date_match is None:
    raise RuntimeError("CITATION.cff has no release date")

RELEASE_DATE = os.environ.get("SEAMS_RELEASE_DATE", _date_match.group(1))
if re.fullmatch(r"\d{4}-\d{2}-\d{2}", RELEASE_DATE) is None:
    raise RuntimeError("SEAMS_RELEASE_DATE must be YYYY-MM-DD")
