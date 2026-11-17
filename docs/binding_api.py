#!/usr/bin/env python3
# Copyright (c) d-SEAMS developers
"""Write the language API pages from the engine headers and the binding repositories.

Each page has the native calls for that language and the C functions those
calls use. ``--check`` fails when a documented name is not in the source,
when a public name is missing from the page, or when a page on disk is not
the page this script would write.
"""

import argparse
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DOCS = ROOT / "docs"
sys.path.insert(0, str(DOCS))
import render_capi  # noqa: E402

OUT = DOCS / "source" / "reference"
HEADER = ROOT / "include" / "seams.h"
LANG_HEADER = ROOT / "bindings" / "cxx" / "seams_lang.h"


def c_names():
    """Every seams_* and SEAMS_* name the C headers export."""
    found = set()
    for path in (HEADER, LANG_HEADER):
        text = path.read_text(encoding="utf-8")
        found.update(re.findall(r"\bseams_[A-Za-z0-9_]+\b", text))
        found.update(re.findall(r"\bSEAMS_[A-Z0-9_]+\b", text))
    found.discard("SEAMS_API")
    return {name for name in found if not name.endswith("_H_")}


def _rst(text):
    """Escape characters that reStructuredText treats as markup."""
    return text.replace("\\", "\\\\").replace("*", "\\*").replace("`", "\\`")


def _heading(title, level):
    char = "=-~"[level]
    return title + "\n" + (char * len(title))


def _fence(lang, body):
    return ".. code-block:: {}\n\n   {}\n".format(lang, body.replace("\n", "\n   "))


def render_page(slug, title, intro, symbols, native_lang):
    """One reference page. ``symbols`` are dicts with name, signature, summary, params, returns, c_api."""
    lines = [_heading(title, 0), "", intro.rstrip(), ""]
    lines.append(_heading("Native API", 1))
    lines.append("")
    lines.append(
        "Generated from the sources named below. "
        "A name on this page is a public call in that language."
    )
    lines.append("")
    for item in symbols:
        lines.append(_heading(item["name"], 2))
        lines.append("")
        lines.append(".. api-symbol: {}".format(item["name"]))
        lines.append("")
        if item.get("signature"):
            lines.append(_fence(native_lang, item["signature"]))
        if item.get("summary"):
            lines.append(_rst(item["summary"]))
            lines.append("")
        if item.get("params"):
            lines.append("Parameters")
            lines.append("")
            for name, detail in item["params"]:
                lines.append("- **{}**: {}".format(_rst(name), _rst(detail or "see the signature")))
            lines.append("")
        if item.get("returns"):
            lines.append("Returns")
            lines.append("")
            lines.append(_rst(item["returns"]))
            lines.append("")
        if item.get("example"):
            lines.append("Example")
            lines.append("")
            lines.append(_fence(native_lang, item["example"]))
    lines.append(_heading("C API", 1))
    lines.append("")
    lines.append(
        "Each native call and the ``seams.h`` or ``seams_lang.h`` function it uses. "
        "The C column is the stable ABI."
    )
    lines.append("")
    lines.append(".. list-table::")
    lines.append("   :header-rows: 1")
    lines.append("   :widths: 40 60")
    lines.append("")
    lines.append("   * - Native")
    lines.append("     - C")
    for item in symbols:
        mapped = item["c_api"]
        if mapped:
            cell = ", ".join("``{}``".format(name) for name in mapped)
        else:
            cell = "no direct C call"
        lines.append("   * - ``{}``".format(item["name"]))
        lines.append("     - {}".format(cell))
    lines.append("")
    return "\n".join(lines)


def c_symbols():
    records = render_capi.parse_header(HEADER.read_text(encoding="utf-8"))
    symbols = []
    seen = set()
    for record in records:
        name = record["name"]
        if name in seen:
            continue
        seen.add(name)
        symbols.append(
            {
                "name": name,
                "signature": record["signature"],
                "summary": record["brief"],
                "params": record["params"],
                "returns": record["returned"],
                "c_api": [name],
            }
        )
    symbols.sort(key=lambda item: item["name"])
    return symbols


def _cpp_calls(body):
    return sorted(set(re.findall(r"\b(seams_[A-Za-z0-9_]+)\s*\(", body)))


def _cpp_brief(text, end):
    window = text[max(0, end - 800) : end]
    matches = list(re.finditer(r"/\*\*(.*?)\*/", window, re.S))
    if not matches:
        return ""
    last = matches[-1]
    if window[last.end() :].strip():
        return ""
    return render_capi._parse_comment(last.group(1))["brief"]


def cpp_symbols():
    """Public C++ types and the methods that call the C API."""
    text = (ROOT / "include" / "seams.hpp").read_text(encoding="utf-8")
    symbols = []
    seen = set()
    for match in re.finditer(r"^(class|struct)\s+(\w+)\b[^{;]*\{", text, re.M):
        name = match.group(2)
        if name == "Handle":
            continue
        start = match.end() - 1
        depth = 0
        body_end = start
        for index in range(start, len(text)):
            char = text[index]
            if char == "{":
                depth += 1
            elif char == "}":
                depth -= 1
                if depth == 0:
                    body_end = index
                    break
        body = text[start : body_end + 1]
        brief = _cpp_brief(text, match.start())
        qual = "seams::" + name
        if qual not in seen:
            seen.add(qual)
            symbols.append(
                {
                    "name": qual,
                    "signature": match.group(0).split("{")[0].strip(),
                    "summary": brief,
                    "params": [],
                    "returns": "",
                    "c_api": _cpp_calls(body),
                }
            )
        access = "private" if match.group(1) == "class" else "public"
        local_calls = {}
        pending = []
        for method in re.finditer(
            r"(public|private|protected)\s*:|(?P<sig>(?:\[\[nodiscard\]\]\s+)?[\w:<>,\s\*&]+?\s+(?P<name>~?\w+)\s*\([^;{}]*\)\s*(?:const)?\s*(?:noexcept)?\s*(?:override)?\s*\{)",
            body,
        ):
            if method.group(1):
                access = method.group(1)
                continue
            method_name = method.group("name")
            if method_name in {"if", "for", "while", "switch"} or method_name.startswith("operator"):
                continue
            chunk_start = method.end() - 1
            depth = 0
            chunk_end = chunk_start
            for index in range(chunk_start, len(body)):
                char = body[index]
                if char == "{":
                    depth += 1
                elif char == "}":
                    depth -= 1
                    if depth == 0:
                        chunk_end = index
                        break
            chunk = body[chunk_start : chunk_end + 1]
            local_calls[method_name] = _cpp_calls(chunk)
            pending.append((access, method_name, method, chunk))
        def resolved(name, stack):
            if name in stack:
                return set()
            found = set(local_calls.get(name, ()))
            chunk = next((item[3] for item in pending if item[1] == name), "")
            for callee in re.findall(r"\b(\w+)\s*\(", chunk):
                if callee in local_calls:
                    found.update(resolved(callee, stack | {name}))
            return found
        for access_name, method_name, method, chunk in pending:
            if access_name != "public":
                continue
            qual_method = qual + "::" + method_name
            if qual_method in seen:
                continue
            seen.add(qual_method)
            signature = " ".join(method.group("sig").split())
            symbols.append(
                {
                    "name": qual_method,
                    "signature": signature[:-1].strip(),
                    "summary": _cpp_brief(body, method.start()),
                    "params": [],
                    "returns": "",
                    "c_api": sorted(resolved(method_name, set())),
                }
            )
    for match in re.finditer(
        r"^(?:\[\[nodiscard\]\] )?inline [\w:<>,\s\*&]+ (\w+)\(([^;]*)\)\s*\{",
        text,
        re.M,
    ):
        name = match.group(1)
        window = text[max(0, match.start() - 400) : match.start()]
        if "@ingroup seams_cpp" not in window:
            continue
        qual = "seams::" + name
        if qual in seen:
            continue
        seen.add(qual)
        body_start = text.find("{", match.end() - 1)
        depth = 0
        body_end = body_start
        for index in range(body_start, len(text)):
            if text[index] == "{":
                depth += 1
            elif text[index] == "}":
                depth -= 1
                if depth == 0:
                    body_end = index
                    break
        symbols.append(
            {
                "name": qual,
                "signature": " ".join(match.group(0).split())[:-1].strip(),
                "summary": _cpp_brief(text, match.start()),
                "params": [],
                "returns": "",
                "c_api": _cpp_calls(text[body_start : body_end + 1]),
            }
        )
    symbols.sort(key=lambda item: item["name"])
    return symbols


LANGUAGES = {
    "c-api": {
        "title": "C",
        "intro": (
            "``include/seams.h`` is the native C API and the ABI the other "
            "languages call. The pages under :doc:`/capi/index` group the same "
            "declarations. This page lists every ``SEAMS_API`` function the "
            "header declares, with the signature and the parameter notes from "
            "that header.\n\n"
            "The C column is the function itself: this is the C API."
        ),
        "lang": "c",
        "load": c_symbols,
    },
    "cpp-api": {
        "title": "C++",
        "intro": (
            "``include/seams.hpp`` is the native C++ API: move-only handles in "
            "namespace ``seams``, with results returned by value. "
            ":doc:`/stable/cpp` is the narrative for that header. "
            "This page lists the public types and methods, and the ``seams.h`` "
            "functions each method calls."
        ),
        "lang": "cpp",
        "load": cpp_symbols,
    },
}


def write_pages(selected):
    OUT.mkdir(parents=True, exist_ok=True)
    written = {}
    for slug in selected:
        spec = LANGUAGES[slug]
        symbols = spec["load"]()
        text = render_page(slug, spec["title"], spec["intro"], symbols, spec["lang"])
        path = OUT / (slug + ".rst")
        path.write_text(text, encoding="utf-8")
        written[slug] = (path, symbols, text)
    return written


def problems(selected, written=None):
    """Missing symbols, names that are not in the source, and stale pages."""
    known = c_names()
    errors = []
    if written is None:
        written = {}
        for slug in selected:
            spec = LANGUAGES[slug]
            symbols = spec["load"]()
            text = render_page(slug, spec["title"], spec["intro"], symbols, spec["lang"])
            written[slug] = (OUT / (slug + ".rst"), symbols, text)
    for slug, (path, symbols, fresh) in written.items():
        names = [item["name"] for item in symbols]
        if len(names) != len(set(names)):
            errors.append("{} has a repeated symbol".format(slug))
        documented = re.findall(r"(?m)^\.\. api-symbol: (\S+)\s*$", fresh)
        if documented != names:
            errors.append("{} symbol block does not match the source list".format(slug))
        for item in symbols:
            for c_name in item["c_api"]:
                if c_name not in known:
                    errors.append(
                        "{} maps {} to unknown C symbol {}".format(slug, item["name"], c_name)
                    )
        if not names:
            errors.append("{} has no public symbols".format(slug))
        if path.is_file():
            have = path.read_text(encoding="utf-8")
            if have != fresh:
                errors.append("{} is stale; regenerate it".format(path.name))
            on_disk = re.findall(r"(?m)^\.\. api-symbol: (\S+)\s*$", have)
            for name in on_disk:
                if name not in names:
                    errors.append("documented symbol does not exist: {}: {}".format(slug, name))
            for name in names:
                if name not in on_disk:
                    errors.append("undocumented {}: {}".format(slug, name))
    return errors


def main(argv=None):
    parser = argparse.ArgumentParser(description="Generate the language API pages")
    parser.add_argument("--check", action="store_true")
    parser.add_argument("languages", nargs="*", default=sorted(LANGUAGES))
    args = parser.parse_args(argv)
    selected = args.languages or sorted(LANGUAGES)
    unknown = [name for name in selected if name not in LANGUAGES]
    if unknown:
        print("unknown language " + " ".join(unknown), file=sys.stderr)
        return 1
    if args.check:
        errors = problems(selected)
        if errors:
            print("\n".join(errors), file=sys.stderr)
            return 1
        count = sum(len(LANGUAGES[slug]["load"]()) for slug in selected)
        print("binding api ok ({} symbols)".format(count))
        return 0
    written = write_pages(selected)
    errors = problems(selected, written)
    if errors:
        print("\n".join(errors), file=sys.stderr)
        return 1
    print("wrote", " ".join(slug for slug in selected))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
