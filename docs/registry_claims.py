#!/usr/bin/env python3
"""Shared install-claim parser and live package-registry lookups."""

from __future__ import annotations

import json
import re
import shlex
import urllib.error
import urllib.parse
import urllib.request

USER_AGENT = "d-SEAMS-reader-text-check/1.0 (+https://docs.dseams.info/)"


class ClaimParseError(ValueError):
    """An install command is present but malformed."""


class RegistryError(RuntimeError):
    """A live registry could not be checked."""


def fetch(url: str, *, missing_ok: bool = False) -> bytes | None:
    request = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
    try:
        with urllib.request.urlopen(request, timeout=20) as response:
            return response.read()
    except urllib.error.HTTPError as error:
        if missing_ok and error.code == 404:
            return None
        raise RegistryError(f"{url}: HTTP {error.code}") from error
    except (OSError, urllib.error.URLError) as error:
        raise RegistryError(f"{url}: {error}") from error


def fetch_json(url: str, *, missing_ok: bool = False):
    payload = fetch(url, missing_ok=missing_ok)
    if payload is None:
        return None
    try:
        return json.loads(payload)
    except (UnicodeDecodeError, json.JSONDecodeError) as error:
        raise RegistryError(f"{url}: invalid JSON: {error}") from error


def live_versions(manager: str, package: str) -> set[str]:
    """Return versions published by the registry used by an install command."""
    quoted = urllib.parse.quote(package)
    if manager == "pip":
        data = fetch_json(f"https://pypi.org/pypi/{quoted}/json", missing_ok=True)
        return set(data.get("releases", ())) if data else set()
    if manager == "cargo":
        data = fetch_json(f"https://crates.io/api/v1/crates/{quoted}", missing_ok=True)
        return {
            item["num"]
            for item in (data or {}).get("versions", ())
            if not item.get("yanked") and item.get("num")
        }
    if manager == "conda":
        data = fetch_json(
            f"https://api.anaconda.org/package/conda-forge/{quoted}",
            missing_ok=True,
        )
        if not data:
            return set()
        versions = set(data.get("versions", ()))
        versions.update(
            item["version"] for item in data.get("files", ()) if item.get("version")
        )
        return versions
    if manager == "luarocks":
        search = fetch(
            "https://luarocks.org/search?q=" + urllib.parse.quote_plus(package)
        ).decode("utf-8", errors="replace")
        module = re.search(
            r'href=["\'](/modules/[^/"\']+/' + re.escape(package) + r')["\']',
            search,
            re.IGNORECASE,
        )
        if not module:
            return set()
        page = fetch("https://luarocks.org" + module.group(1)).decode(
            "utf-8", errors="replace"
        )
        return set(re.findall(re.escape(module.group(1)) + r'/([^/"\']+)["\']', page))
    if manager == "julia":
        initial = urllib.parse.quote(package[:1].upper())
        data = fetch(
            "https://raw.githubusercontent.com/JuliaRegistries/General/master/"
            f"{initial}/{quoted}/Versions.toml",
            missing_ok=True,
        )
        if data is None:
            return set()
        return set(re.findall(r'^\["([^"]+)"\]', data.decode(), re.MULTILINE))
    if manager == "r":
        versions = set()
        description = fetch(
            f"https://cran.r-project.org/web/packages/{quoted}/DESCRIPTION",
            missing_ok=True,
        )
        if description:
            match = re.search(r"^Version:\s*(\S+)", description.decode(), re.MULTILINE)
            if match:
                versions.add(match.group(1))
        data = fetch_json(
            f"https://d-seams.r-universe.dev/api/packages/{quoted}",
            missing_ok=True,
        )
        if data and data.get("Version"):
            versions.add(data["Version"])
        return versions
    if manager == "fpm":
        data = fetch_json(
            "https://raw.githubusercontent.com/fortran-lang/"
            "fpm-registry/master/index.json"
        )
        records = data.get("packages", {}).get(package.lower(), {})
        versions = {key for key in records if key != "latest"}
        latest = records.get("latest", {})
        if latest.get("version"):
            versions.add(latest["version"])
        return versions
    raise RegistryError(f"unsupported package manager: {manager}")


def claims(line: str):
    """Yield (manager, package, version) install claims found in one line."""
    plain = line.replace("~", "").replace("`", "")

    for match in re.finditer(
        r"\bpip(?:3)?\s+install(?:\s+-\S+)*\s+['\"]?([A-Za-z0-9_.-]+)"
        r"(?:==([A-Za-z0-9_.+-]+))?",
        plain,
    ):
        package = match.group(1).rstrip(".,;:")
        if package:
            yield "pip", package, (match.group(2) or "*").rstrip(".,;:")

    match = re.search(r"\b(?:conda|mamba)\s+install\b(.*)", plain)
    if match:
        try:
            tokens = shlex.split(match.group(1))
        except ValueError as error:
            raise ClaimParseError(
                f"malformed conda install command: {match.group(0)} ({error})"
            ) from error
        skip = False
        for token in tokens:
            if skip:
                skip = False
                continue
            if token in {"-c", "--channel"}:
                skip = True
                continue
            if token.startswith("-"):
                continue
            package, sep, version = token.partition("=")
            yield "conda", package.rstrip(".,;:"), (version if sep else "*").rstrip(".,;:")

    for match in re.finditer(
        r"\bluarocks\s+install\s+([A-Za-z0-9_.-]+)(?:\s+([A-Za-z0-9_.+-]+))?",
        plain,
    ):
        yield "luarocks", match.group(1).rstrip(".,;:"), (match.group(2) or "*").rstrip(".,;:")

    if not re.search(r"\bcargo\s+add\b[^\n]*\s--git(?:\s|=)", plain):
        for match in re.finditer(
            r"\bcargo\s+add(?:\s+--\S+(?:\s+\S+)?)?\s+([A-Za-z0-9_.-]+)"
            r"(?:@([A-Za-z0-9_.+-]+))?",
            plain,
        ):
            yield "cargo", match.group(1).rstrip(".,;:"), (match.group(2) or "*").rstrip(".,;:")

    for call in re.finditer(r"\bPkg\.add\(([^)]*)\)", plain):
        body = call.group(1)
        if "url=" in body:
            continue
        name = re.search(r'name\s*=\s*["\']([^"\']+)', body)
        version = re.search(r'version\s*=\s*["\']([^"\']+)', body)
        if name:
            yield "julia", name.group(1), version.group(1) if version else "*"
        else:
            short = re.fullmatch(r'\s*["\']([^"\']+)["\']\s*', body)
            if short:
                yield "julia", short.group(1), "*"

    for match in re.finditer(r"\binstall\.packages\(\s*[\"']([^\"']+)[\"']", plain):
        yield "r", match.group(1), "*"

    for match in re.finditer(
        r"\bfpm\s+(?:install|package)\s+[`'\"]?([A-Za-z0-9_.-]+)"
        r"(?:@([A-Za-z0-9_.+-]+))?",
        plain,
        re.IGNORECASE,
    ):
        yield "fpm", match.group(1).rstrip(".,;:"), (match.group(2) or "*").rstrip(".,;:")
