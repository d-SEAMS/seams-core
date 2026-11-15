#!/usr/bin/env python3
# Copyright (c) d-SEAMS developers
"""Every exported C function is documented in seams.h and on the C API page.

Usage: capi_docs.py [--write] <seams.h> <capi_symbols.txt> <docs/api/c.md>

--write first rewrites the page's section of every listed function from
the header, then checks.

In seams.h each function listed in capi_symbols.txt needs a doc comment
directly above its declaration with @brief, an @param for every parameter
and none for a parameter it does not have, @return unless it returns void,
@par Thread safety, and @par Ownership when a pointer goes in or comes out.
On the page, the symbol reference block needs a section for the function
that carries the same purpose, parameters, return, ownership and thread
safety text, as render_c_api.py writes it.
"""

import re
import sys
from dataclasses import dataclass, field
from pathlib import Path

BEGIN = "<!-- API_SYMBOLS_BEGIN -->"
END = "<!-- API_SYMBOLS_END -->"
DECL = re.compile(
    r"(?:(/\*\*(?:(?!\*/).)*?\*/)\s*)?"
    r"(?:SEAMS_DEPRECATED\((?:[^()]|\([^()]*\))*\)\s*)?"
    r"SEAMS_API\s+([^;(]*?)\b(seams_[A-Za-z0-9_]+)\s*\(([^;]*?)\)\s*;",
    re.S,
)


@dataclass
class Function:
    name: str
    returns_type: str
    declaration: str
    parameters: list
    documented: bool = False
    brief: str = ""
    details: str = ""
    params: dict = field(default_factory=dict)
    returns: str = ""
    sections: dict = field(default_factory=dict)

    @property
    def takes_pointer(self):
        return "*" in self.returns_type or any("*" in p for p in self.parameters)

    @property
    def returns_void(self):
        return self.returns_type.split() == ["void"]

    def param_names(self):
        names = []
        for decl in self.parameters:
            found = re.findall(r"([A-Za-z_][A-Za-z0-9_]*)\s*(?:\[[^\]]*\])?\s*$", decl)
            names.append(found[0] if found else decl)
        return names


def _comment_lines(comment):
    body = comment[3:-2]
    lines = []
    for row in body.splitlines():
        row = row.strip()
        if row.startswith("*"):
            row = row[1:]
            if row.startswith(" "):
                row = row[1:]
        lines.append(row.rstrip())
    return lines


def _clean(text):
    text = re.sub(r"@p\s+(\w+)", r"`\1`", text)
    text = re.sub(r"::(SEAMS_\w+)", r"\1", text)
    return " ".join(text.split())


def _parse_comment(fn, comment):
    fn.documented = True
    current = None
    buf = []
    code = []
    in_code = False

    def flush():
        text = _clean(" ".join(buf))
        if current is None:
            if text:
                fn.details = (fn.details + " " + text).strip()
        elif current[0] == "brief":
            fn.brief = text
        elif current[0] == "param":
            fn.params[current[1]] = text
        elif current[0] == "return":
            fn.returns = text
        elif current[0] == "par":
            fn.sections[current[1]] = text
        buf.clear()

    for row in _comment_lines(comment):
        if in_code:
            if row.strip() == "@endcode":
                in_code = False
                fn.sections["Example code"] = "\n".join(code)
            else:
                code.append(row)
            continue
        if row.strip() == "@code":
            in_code = True
            code = []
            continue
        tag = re.match(r"@(brief|param|return|retval|par)\b(\[[a-z,]+\])?\s*(.*)$", row.strip())
        if tag:
            flush()
            kind, _, rest = tag.groups()
            if kind == "brief":
                current = ("brief",)
                buf.append(rest)
            elif kind == "param":
                pname, _, text = rest.partition(" ")
                current = ("param", pname)
                buf.append(text)
            elif kind in ("return", "retval"):
                current = ("return",)
                buf.append(rest)
            else:
                current = ("par", rest.strip())
            continue
        if not row.strip() and current is not None and current[0] in ("brief",):
            flush()
            current = None
            continue
        buf.append(row)
    flush()


def functions(header_text):
    """Every SEAMS_API function in the header, with its doc comment parsed."""
    text = header_text
    found = {}
    for match in DECL.finditer(text):
        comment, ret, name, params = match.groups()
        decl = " ".join(("SEAMS_API " + ret + name + "(" + params + ");").split())
        plist = [" ".join(p.split()) for p in params.split(",")]
        plist = [p for p in plist if p and p != "void"]
        fn = Function(name, " ".join(ret.split()), decl.replace("SEAMS_API ", "", 1), plist)
        # DECL only takes a comment that ends right above the declaration.
        if comment is not None:
            _parse_comment(fn, comment)
        found[name] = fn
    return found


def header_problems(funcs, names):
    errors = []
    for name in names:
        fn = funcs.get(name)
        if fn is None:
            errors.append("{}: not declared in seams.h".format(name))
            continue
        if not fn.documented:
            errors.append("{}: no doc comment above the declaration".format(name))
            continue
        if not fn.brief:
            errors.append("{}: no @brief".format(name))
        wanted = fn.param_names()
        for pname in wanted:
            if not fn.params.get(pname):
                errors.append("{}: no @param for {}".format(name, pname))
        for pname in fn.params:
            if pname not in wanted:
                errors.append("{}: @param {} names no parameter".format(name, pname))
        if not fn.returns_void and not fn.returns:
            errors.append("{}: no @return".format(name))
        if not fn.sections.get("Thread safety"):
            errors.append("{}: no @par Thread safety".format(name))
        if fn.takes_pointer and not fn.sections.get("Ownership"):
            errors.append("{}: takes or returns a pointer and has no @par Ownership".format(name))
    return errors


def page_sections(page_text):
    if BEGIN not in page_text or END not in page_text:
        return None
    body = page_text.split(BEGIN, 1)[1].split(END, 1)[0]
    found = {}
    for chunk in re.split(r"^### ", body, flags=re.M)[1:]:
        title, _, rest = chunk.partition("\n")
        found[title.strip()] = rest
    return found


def page_problems(funcs, names, page_text):
    sections = page_sections(page_text)
    if sections is None:
        return ["the C API page has no symbol reference block"]
    errors = []
    for name in names:
        fn = funcs.get(name)
        block = sections.get(name)
        if block is None:
            errors.append("{}: no section on the C API page".format(name))
            continue
        if fn is None or not fn.documented:
            continue
        flat = " ".join(block.split())
        wanted = ["Purpose: " + fn.brief]
        wanted += ["- `{}`: {}".format(p, fn.params.get(p, "")) for p in fn.param_names()]
        if not fn.returns_void:
            wanted.append("Returns: " + fn.returns)
        if fn.sections.get("Ownership"):
            wanted.append("Ownership: " + fn.sections["Ownership"])
        wanted.append("Thread safety: " + fn.sections.get("Thread safety", ""))
        for item in wanted:
            if " ".join(item.split()) not in flat:
                label = item.split(":", 1)[0]
                errors.append("{}: the C API page lacks {} as documented in seams.h".format(name, label))
    return errors


def render_function(fn):
    """Symbol reference section of one function, as markdown chunks."""
    out = ["### {}\n".format(fn.name), "Purpose: {}\n".format(fn.brief)]
    if fn.details:
        out.append(fn.details + "\n")
    out.append("```c\n{}\n```\n".format(fn.declaration))
    names = fn.param_names()
    if names:
        out.append("Parameters: {}.\n".format(", ".join("`{}`".format(p) for p in names)))
        out.append("\n".join("- `{}`: {}".format(p, fn.params.get(p, "")) for p in names) + "\n")
    else:
        out.append("Parameters: none.\n")
    out.append("Returns: {}\n".format("nothing." if fn.returns_void else fn.returns))
    for title, text in fn.sections.items():
        if title not in ("Example", "Example code") and text:
            out.append("{}: {}\n".format(title, text))
    out.append("Units: as given for each parameter above.\n")
    if "Example code" in fn.sections:
        out.append("Example:\n")
        out.append("```c\n{}\n```\n".format(fn.sections["Example code"]))
    out.append("The tour call is in `docs/examples/c_api_tour.c`.\n")
    return out


def rewrite_page(page_text, funcs, names):
    """Replace the sections of the listed functions inside the symbol block."""
    head, rest = page_text.split(BEGIN, 1)
    body, tail = rest.split(END, 1)
    parts = re.split(r"^(?=### )", body, flags=re.M)
    seen = set()
    out = []
    for part in parts:
        title = part[4:].partition("\n")[0].strip() if part.startswith("### ") else None
        fn = funcs.get(title) if title in names else None
        if fn is not None and fn.documented:
            out.append("\n".join(render_function(fn)) + "\n")
            seen.add(title)
        else:
            out.append(part)
    for name in names:
        fn = funcs.get(name)
        if name not in seen and fn is not None and fn.documented:
            out.append("\n".join(render_function(fn)) + "\n")
    return head + BEGIN + "".join(out) + END + tail


def listed_names(path):
    rows = Path(path).read_text().splitlines()
    return [row.strip() for row in rows if row.strip() and not row.startswith("#")]


def problems(header_text, names, page_text):
    funcs = functions(header_text)
    return header_problems(funcs, names) + page_problems(funcs, names, page_text)


def main(argv):
    write = len(argv) > 1 and argv[1] == "--write"
    if write:
        argv = argv[:1] + argv[2:]
    if len(argv) != 4:
        print(__doc__, file=sys.stderr)
        return 2
    names = listed_names(argv[2])
    if write:
        page = Path(argv[3])
        page.write_text(rewrite_page(page.read_text(), functions(Path(argv[1]).read_text()), names))
    errors = problems(Path(argv[1]).read_text(), names, Path(argv[3]).read_text())
    if errors:
        print("\n".join(errors), file=sys.stderr)
        return 1
    print("capi docs complete ({} functions)".format(len(names)))
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
