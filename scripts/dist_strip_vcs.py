# Copyright (c) d-SEAMS developers
"""Drop VCS metadata from subprojects in a `meson dist` tree.

A wrap-git checkout whose meson.build comes from a packagefile is not
tracked in its own repository, so `meson dist` copies the directory,
.git included, instead of archiving it.
"""

import os
import shutil
from pathlib import Path

root = Path(os.environ["MESON_DIST_ROOT"]) / "subprojects"
if root.is_dir():
    for child in root.iterdir():
        vcs = child / ".git"
        if vcs.is_dir():
            shutil.rmtree(vcs)
        elif vcs.exists():
            vcs.unlink()
