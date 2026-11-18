#!/usr/bin/env python3
"""Validate documented install claims against live package registries."""

from __future__ import annotations

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from check_reader_text import source_pages
from registry_claims import ClaimParseError, RegistryError, claims, live_versions


def check_registry_pages(pages, resolve_versions=live_versions):
    errors = []
    registry_cache = {}
    registry_claims = 0
    scanned_files = 0

    for path, rel in pages:
        scanned_files += 1
        text = path.read_text(encoding="utf-8", errors="replace")
        for number, line in enumerate(text.splitlines(), 1):
            try:
                line_claims = list(claims(line))
            except ClaimParseError as error:
                errors.append(f"{rel}:{number}: malformed install claim: {error}")
                continue
            for manager, package, version in line_claims:
                registry_claims += 1
                key = (manager, package.lower())
                if key not in registry_cache:
                    try:
                        registry_cache[key] = resolve_versions(*key)
                    except RegistryError as error:
                        registry_cache[key] = error
                available = registry_cache[key]
                if isinstance(available, RegistryError):
                    errors.append(
                        f"{rel}:{number}: live registry check failed for "
                        f"{manager} {package}: {available}"
                    )
                elif not available or (version != "*" and version not in available):
                    errors.append(
                        f"{rel}:{number}: unpublished install claim: "
                        f"{manager} {package} {version}"
                    )
    return errors, scanned_files, registry_claims, len(registry_cache)


def main() -> int:
    errors, scanned_files, registry_claims, package_count = check_registry_pages(
        source_pages()
    )
    if errors:
        print("\n".join(errors), file=sys.stderr)
        return 1
    print(
        f"registry claims ok ({scanned_files} exported sources; "
        f"{registry_claims} registry claims checked live across {package_count} packages)"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
