#!/usr/bin/env python3
# Copyright (c) d-SEAMS developers
"""Fail when a language API page names a symbol the source does not have.

The pages are written by docs/binding_api.py. This check reruns that
extraction. A documented name with no source declaration is stale. A
public declaration with no page entry is missing.
"""

import binding_api


if __name__ == "__main__":
    raise SystemExit(binding_api.main(["--check"]))
