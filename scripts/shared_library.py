# Copyright (c) d-SEAMS developers
"""Locate the shared library Meson built, without assuming one suffix.

``shared_library.full_path()`` is passed in an environment variable when
the tests are launched from Meson. That file wins. Otherwise the
candidate that exists is used, and ``lib<stem>.so`` is tried before
``lib<stem>.dylib``, so a Linux build still resolves to the same path.
"""

import os
from pathlib import Path


def library_names(stem):
    names = ["lib" + stem + ".so", "lib" + stem + ".dylib"]
    if stem == "seams":
        names.append("seams.dll")
    else:
        names.append("lib" + stem + ".dll")
    return names


def built_library(directory, stem, env_var=None):
    """Return the built library, or None when it is not there."""
    if env_var:
        named = os.environ.get(env_var, "").strip()
        if named and Path(named).is_file():
            return Path(named)
    directory = Path(directory)
    for name in library_names(stem):
        candidate = directory / name
        if candidate.is_file():
            return candidate
    return None


def require_library(directory, stem, env_var=None):
    found = built_library(directory, stem, env_var)
    if found is None:
        raise SystemExit("lib" + stem + " was not built in " + str(directory))
    return found
