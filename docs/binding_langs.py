# Copyright (c) d-SEAMS developers
"""Parsers for the language bindings and the seams command.

The clones live under docs/.cache/bindings and are refreshed from
github.com/d-SEAMS when this module loads a language.
"""

import ast
import os
import re
import shutil
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
CACHE = ROOT / "docs" / ".cache" / "bindings"

REPOS = {
    "PydSEAMSlib": "https://github.com/d-SEAMS/PydSEAMSlib.git",
    "Dseams.jl": "https://github.com/d-SEAMS/Dseams.jl.git",
    "dseams-r": "https://github.com/d-SEAMS/dseams-r.git",
    "dseams-rs": "https://github.com/d-SEAMS/dseams-rs.git",
    "dseams-fortran": "https://github.com/d-SEAMS/dseams-fortran.git",
    "yodaStruct": "https://github.com/d-SEAMS/yodaStruct.git",
}


def checkout(name):
    dest = CACHE / name
    dest.parent.mkdir(parents=True, exist_ok=True)
    # The offline tarball has no network. meson dist vendors these clones
    # into docs/.cache/bindings and strips .git, so a fetch cannot refresh them.
    if os.environ.get("SEAMS_NO_NETWORK") == "1":
        if dest.is_dir() and any(dest.iterdir()):
            return dest
        raise RuntimeError(
            "{} is not in docs/.cache/bindings and SEAMS_NO_NETWORK=1 forbids cloning".format(name)
        )
    if not (dest / ".git").is_dir():
        # A dist tree keeps the sources and drops .git. An online run replaces
        # that snapshot so the check still tracks the binding heads.
        if dest.exists():
            shutil.rmtree(dest)
        subprocess.check_call(["git", "clone", "--quiet", "--depth", "1", REPOS[name], str(dest)])
    else:
        subprocess.check_call(["git", "-C", str(dest), "fetch", "--quiet", "--depth", "1", "origin", "HEAD"])
        subprocess.check_call(["git", "-C", str(dest), "reset", "--quiet", "--hard", "FETCH_HEAD"])
    return dest


def _calls(text):
    return sorted(set(re.findall(r"\b(seams_[A-Za-z0-9_]+)\s*\(", text)))


def _numpy_doc(doc):
    if not doc:
        return "", [], "", ""
    lines = doc.strip().splitlines()
    summary = []
    params = []
    returns = []
    example = []
    section = "summary"
    current = None
    for line in lines:
        stripped = line.strip()
        if stripped in {"Parameters", "Returns", "Examples", "Example"} and set(line[len(line) - len(stripped) - 1 :].strip()) <= set("-"):
            section = stripped.lower()
            current = None
            continue
        if section == "summary" and stripped and not stripped.startswith("---"):
            summary.append(stripped)
        elif section == "parameters":
            if stripped and not line.startswith(" ") and ":" in stripped and not stripped.startswith("---"):
                name, _, detail = stripped.partition(":")
                current = [name.strip(), detail.strip()]
                params.append(current)
            elif current is not None and stripped:
                current[1] = (current[1] + " " + stripped).strip()
        elif section.startswith("return"):
            if stripped and not stripped.startswith("---"):
                returns.append(stripped)
        elif section.startswith("example"):
            example.append(line[4:] if line.startswith("    ") else line)
    return " ".join(summary), params, " ".join(returns), "\n".join(example).strip()


def _signature(node):
    args = []
    positional = node.args.args
    defaults = [None] * (len(positional) - len(node.args.defaults)) + list(node.args.defaults)
    for arg, default in zip(positional, defaults):
        if arg.arg in {"self", "cls"}:
            continue
        if default is None:
            args.append(arg.arg)
        else:
            try:
                rendered = ast.unparse(default)
            except Exception:
                rendered = "..."
            args.append("{}={}".format(arg.arg, rendered))
    if node.args.vararg:
        args.append("*" + node.args.vararg.arg)
    if node.args.kwarg:
        args.append("**" + node.args.kwarg.arg)
    return args


def _py_item(qual, node, source, core):
    doc = ast.get_docstring(node) or ""
    summary, params, returns, example = _numpy_doc(doc)
    calls = set(core.get(node.name, ()))
    segment = ast.get_source_segment(source, node) or ""
    calls.update(_calls(segment))
    for callee in re.findall(r"\.(\w+)\s*\(", segment):
        calls.update(core.get(callee, ()))
    args = _signature(node) if isinstance(node, ast.FunctionDef) else []
    signature = "{}({})".format(qual, ", ".join(args)) if args or isinstance(node, ast.FunctionDef) else qual
    return {
        "name": qual,
        "signature": signature,
        "summary": summary,
        "params": params,
        "returns": returns,
        "example": example,
        "c_api": sorted(calls),
    }


def _core_map(text):
    parts = re.split(r'\.def\(\s*"(\w+)"', text)
    found = {}
    for index in range(1, len(parts), 2):
        found.setdefault(parts[index], set()).update(_calls(parts[index + 1][:2000]))
    for match in re.finditer(r'm\.def\(\s*"(\w+)"\s*,', text):
        found.setdefault(match.group(1), set()).update(_calls(text[match.end() : match.end() + 500]))
    return found


def python_symbols():
    root = checkout("PydSEAMSlib")
    core = _core_map((root / "src" / "dseams" / "_core.cpp").read_text(encoding="utf-8"))
    symbols = []
    seen = set()
    for path in sorted((root / "src").rglob("*.py")):
        if any(part.startswith("_") and part != "_core.cpp" for part in path.parts):
            if path.name.startswith("_") and path.name != "__init__.py":
                continue
        if path.name.startswith("_"):
            continue
        rel = path.relative_to(root / "src")
        module = ".".join(rel.with_suffix("").parts)
        if module.endswith(".__init__"):
            module = module[: -len(".__init__")]
        source = path.read_text(encoding="utf-8")
        tree = ast.parse(source)
        for node in tree.body:
            if isinstance(node, ast.FunctionDef) and not node.name.startswith("_"):
                qual = module + "." + node.name
                if qual in seen:
                    continue
                seen.add(qual)
                symbols.append(_py_item(qual, node, source, core))
            elif isinstance(node, ast.ClassDef) and not node.name.startswith("_"):
                qual = module + "." + node.name
                if qual not in seen:
                    seen.add(qual)
                    symbols.append(_py_item(qual, node, source, core))
                for child in node.body:
                    if isinstance(child, ast.FunctionDef) and not child.name.startswith("_"):
                        method = qual + "." + child.name
                        if method in seen:
                            continue
                        seen.add(method)
                        symbols.append(_py_item(method, child, source, core))
    symbols.sort(key=lambda item: item["name"])
    return symbols


def _julia_items(text, prefix):
    symbols = []
    pattern = re.compile(r'"""(.*?)"""\s*(?:function\s+)?([A-Za-z][\w!]*)', re.S)
    for match in pattern.finditer(text):
        name = match.group(2)
        if name.startswith("_"):
            continue
        doc = " ".join(match.group(1).split())
        start = match.end()
        stop = text.find("\nend", start)
        body = text[start : stop if stop > 0 else start + 800]
        signature = doc.split(". ", 1)[0]
        if not signature.startswith(name) and "(" not in signature:
            signature = name
        symbols.append(
            {
                "name": prefix + name,
                "signature": signature[:240],
                "summary": doc[:500],
                "params": [],
                "returns": "",
                "example": "",
                "c_api": _calls(body) + sorted(set(re.findall(r":(seams_[A-Za-z0-9_]+)", body))),
            }
        )
    # unique, keep first
    found = {}
    for item in symbols:
        item["c_api"] = sorted(set(item["c_api"]))
        found.setdefault(item["name"], item)
    return sorted(found.values(), key=lambda item: item["name"])


def julia_symbols():
    root = checkout("Dseams.jl")
    items = _julia_items((root / "src" / "Dseams.jl").read_text(encoding="utf-8"), "Dseams.")
    seams = root / "lib" / "Seams" / "src" / "Seams.jl"
    if seams.is_file():
        items.extend(_julia_items(seams.read_text(encoding="utf-8"), "Seams."))
    items.sort(key=lambda item: item["name"])
    return items


def _r_cpp_map(root):
    text = (root / "src" / "dseams.cpp").read_text(encoding="utf-8")
    found = {}
    for match in re.finditer(
        r"\[\[cpp11::register\]\]\s*[\w:<>,\s\*&]+\s+(\w+)\s*\((.*?)\)\s*\{",
        text,
        re.S,
    ):
        name = match.group(1)
        start = match.end() - 1
        depth = 0
        end = start
        for index in range(start, len(text)):
            if text[index] == "{":
                depth += 1
            elif text[index] == "}":
                depth -= 1
                if depth == 0:
                    end = index
                    break
        found[name] = _calls(text[start : end + 1])
    return found


def r_symbols():
    root = checkout("dseams-r")
    cpp_map = _r_cpp_map(root)
    symbols = []
    seen = set()
    for path in sorted((root / "R").glob("*.R")):
        text = path.read_text(encoding="utf-8")
        found = re.findall(
            r"((?:#'[^\n]*\n)+)([A-Za-z.][\w.]*)\s*<-\s*function\s*\(([^)]*)\)(.*?)(?=\n#'|\n[A-Za-z.][\w.]*\s*<-|\Z)",
            text,
            re.S,
        )
        for doc, name, args, body in found:
            if "@export" not in doc or name in seen:
                continue
            seen.add(name)
            title = ""
            params = []
            returns = ""
            example = []
            section = ""
            for line in doc.splitlines():
                if line.startswith("#'"):
                    prose = line[2:].strip()
                    if prose.startswith("@param"):
                        rest = prose[len("@param") :].strip()
                        pname, _, detail = rest.partition(" ")
                        params.append((pname, detail))
                    elif prose.startswith("@return"):
                        returns = prose[len("@return") :].strip()
                    elif prose.startswith("@examples"):
                        section = "ex"
                    elif prose.startswith("@"):
                        section = ""
                    elif section == "ex":
                        example.append(prose)
                    elif not title and prose:
                        title = prose
            calls = set(_calls(body))
            for callee in re.findall(r"\b([A-Za-z][\w]*?)_\s*\(", body):
                calls.update(cpp_map.get(callee + "_", ()))
                calls.update(cpp_map.get(callee, ()))
            symbols.append(
                {
                    "name": name,
                    "signature": "{}({})".format(name, " ".join(args.split())),
                    "summary": title,
                    "params": params,
                    "returns": returns,
                    "example": "\n".join(example).strip(),
                    "c_api": sorted(calls),
                }
            )
    symbols.sort(key=lambda item: item["name"])
    return symbols


def rust_symbols():
    root = checkout("dseams-rs")
    text = (root / "dseams" / "src" / "lib.rs").read_text(encoding="utf-8")
    symbols = []
    seen = set()
    for match in re.finditer(r"(?m)^((?:///.*\n)*)pub (struct|enum|fn) (\w+)", text):
        name = match.group(3)
        if name.startswith("_") or name in seen:
            continue
        seen.add(name)
        docs = " ".join(line[3:].strip() for line in match.group(1).splitlines() if line.startswith("///"))
        kind = match.group(2)
        if kind == "fn":
            end = text.find("{", match.end())
            signature = "pub fn " + name + " " + " ".join(text[match.end() : end].split())
            body_end = end
            depth = 0
            for index in range(end, len(text)):
                if text[index] == "{":
                    depth += 1
                elif text[index] == "}":
                    depth -= 1
                    if depth == 0:
                        body_end = index
                        break
            body = text[end : body_end + 1]
        else:
            signature = "pub {} {}".format(kind, name)
            body = text[match.end() : match.end() + 400]
        calls = sorted(set(re.findall(r"\b((?:ffi|seams_sys)::[A-Za-z0-9_]+|seams_[A-Za-z0-9_]+)\b", body)))
        calls = [name for name in calls if name.startswith("seams_") or name.startswith("ffi::seams") or "seams_" in name]
        plain = []
        for call in calls:
            token = call.split("::")[-1]
            if token.startswith("seams_") or token.startswith("SEAMS_"):
                plain.append(token)
        symbols.append(
            {
                "name": "seams::" + name,
                "signature": signature[:300],
                "summary": docs[:500],
                "params": [],
                "returns": "",
                "example": "",
                "c_api": sorted(set(plain)),
            }
        )
    symbols.sort(key=lambda item: item["name"])
    return symbols


def fortran_symbols():
    root = checkout("dseams-fortran")
    text = (root / "src" / "seams.f90").read_text(encoding="utf-8")
    methods = (root / "src" / "seams_methods.f90").read_text(encoding="utf-8")
    c_text = (root / "src" / "seams_c.f90").read_text(encoding="utf-8")
    bound = set(re.findall(r'bind\(\s*C\s*,\s*name\s*=\s*"([^"]+)"\s*\)', c_text))
    symbols = []
    seen = set()
    for match in re.finditer(r"(?m)^[ \t]*((?:!>[^\n]*\n)+)[ \t]*(?:type,\s*public\s*::\s*(\w+)|(\w+)\s*\()", text):
        docs = " ".join(
            re.sub(r"^\s*!>?\s?", "", line).strip() for line in match.group(1).splitlines()
        )
        name = match.group(2) or match.group(3)
        if not name or name in seen or name in {"if", "do"}:
            continue
        seen.add(name)
        calls = [token for token in bound if token.endswith(name) or name in token]
        if not calls:
            start = match.end()
            calls = _calls(text[start : start + 500])
        symbols.append(
            {
                "name": "seams::" + name,
                "signature": name,
                "summary": docs[:400],
                "params": [],
                "returns": "",
                "example": "",
                "c_api": sorted(set(calls)),
            }
        )
    for match in re.finditer(r"(?m)^[ \t]*procedure,\s*public\s*::\s*(\w+)\s*=>\s*(\w+)", text):
        public, impl = match.group(1), match.group(2)
        qual = "seams::" + public
        if qual in seen:
            continue
        seen.add(qual)
        impl_at = methods.find("module procedure " + impl)
        if impl_at < 0:
            impl_at = methods.find("function " + impl)
        if impl_at < 0:
            impl_at = methods.find("subroutine " + impl)
        chunk = ""
        if impl_at >= 0:
            stop = methods.find("\n  end procedure", impl_at)
            if stop < 0:
                stop = methods.find("\n  end function", impl_at)
            if stop < 0:
                stop = methods.find("\n  end subroutine", impl_at)
            chunk = methods[impl_at : stop if stop > 0 else impl_at + 800]
        calls = set(_calls(chunk))
        for callee in re.findall(r"\bcall\s+(\w+)", chunk):
            at = methods.find("subroutine " + callee)
            if at >= 0:
                stop = methods.find("\n  end subroutine", at)
                calls.update(_calls(methods[at : stop if stop > 0 else at + 800]))
        symbols.append(
            {
                "name": qual,
                "signature": public + " => " + impl,
                "summary": "",
                "params": [],
                "returns": "",
                "example": "",
                "c_api": sorted(calls),
            }
        )
    symbols.sort(key=lambda item: item["name"])
    return symbols


def lua_symbols():
    root = checkout("yodaStruct")
    text = (root / "src" / "lua_api.cpp").read_text(encoding="utf-8")
    header = (root / "src" / "include" / "internal" / "lua_api.hpp").read_text(encoding="utf-8")
    symbols = []
    seen = set()
    for match in re.finditer(r'set_function\(\s*"([^"]+)"', text):
        name = match.group(1)
        if name in seen:
            continue
        seen.add(name)
        window = text[max(0, match.start() - 200) : match.end() + 400]
        calls = _calls(window)
        target = re.search(r'set_function\(\s*"[^"]+"\s*,\s*([A-Za-z0-9_:]+)', match.group(0) + text[match.end() : match.end() + 80])
        if target:
            fn = target.group(1).split("::")[-1]
            at = text.find(fn + "(")
            if at >= 0:
                calls = sorted(set(calls) | set(_calls(text[at : at + 900])))
        brief = ""
        comment = re.search(r"/\*\*(.*?)\*/\s*$", header[max(0, header.find(name) - 400) : header.find(name) + 1], re.S)
        if comment:
            brief = " ".join(re.sub(r"^\s*\*\s?", "", comment.group(1), flags=re.M).split())
        symbols.append(
            {
                "name": "dseams." + name,
                "signature": name,
                "summary": brief[:400],
                "params": [],
                "returns": "",
                "example": "",
                "c_api": calls,
            }
        )
    symbols.sort(key=lambda item: item["name"])
    return symbols


def cli_symbols():
    text = (ROOT / "src" / "seams_cli.cpp").read_text(encoding="utf-8")
    symbols = []
    seen = set()
    for match in re.finditer(r'Option\("(-{1,2}[\w-]+)"', text):
        flag = match.group(1)
        if flag in seen:
            continue
        seen.add(flag)
        help_match = re.search(r'\.help\("([^"]*)"', text[match.end() : match.end() + 400])
        symbols.append(
            {
                "name": "seams " + flag,
                "signature": flag,
                "summary": help_match.group(1) if help_match else "",
                "params": [],
                "returns": "",
                "example": "",
                "c_api": [],
            }
        )
    commands = sorted(set(re.findall(r'cmd == "([a-z0-9][a-z0-9_+-]*)"', text)))
    for command in commands:
        if command in seen:
            continue
        seen.add(command)
        parts = []
        for block in re.finditer(r'cmd == "' + re.escape(command) + r'"', text):
            parts.append(text[block.start() : block.start() + 700])
        symbols.append(
            {
                "name": "seams " + command,
                "signature": "seams " + command + " FILE",
                "summary": "Subcommand {}.".format(command),
                "params": [],
                "returns": "",
                "example": "",
                "c_api": _calls("\n".join(parts)),
            }
        )
    symbols.sort(key=lambda item: item["name"])
    return symbols


def languages():
    return {
        "python-api": {
            "title": "Python",
            "intro": (
                "``pydseamslib`` provides ``import dseams`` and ``import pydseams``. "
                "Signatures, parameters, returns, and examples are read from the "
                "docstrings in `d-SEAMS/PydSEAMSlib <https://github.com/d-SEAMS/PydSEAMSlib>`_. "
                "The C column is taken from the nanobind definitions in ``src/dseams/_core.cpp`` "
                "and from ``seams_*`` calls in the function."
            ),
            "lang": "python",
            "load": python_symbols,
        },
        "julia-api": {
            "title": "Julia",
            "intro": (
                "Documenter-style docstrings from `Dseams.jl <https://github.com/d-SEAMS/Dseams.jl>`_ "
                "and the ``Seams`` package in that repository. "
                "``Dseams`` is the language front. ``Seams`` calls ``libseams`` with ``ccall``."
            ),
            "lang": "julia",
            "load": julia_symbols,
        },
        "r-api": {
            "title": "R",
            "intro": (
                "Rd-style entries for the exported functions in "
                "`dseams-r <https://github.com/d-SEAMS/dseams-r>`_. "
                "The title, ``@param``, ``@return``, and ``@examples`` lines are the "
                "roxygen block above each function."
            ),
            "lang": "r",
            "load": r_symbols,
        },
        "rust-api": {
            "title": "Rust",
            "intro": (
                "Item docs from the ``dseams`` crate, taken from "
                "`dseams/src/lib.rs <https://github.com/d-SEAMS/dseams-rs/blob/main/dseams/src/lib.rs>`_. "
                "``dseams`` is the safe API. ``dseams-sys`` is the raw FFI, re-exported as ``seams::ffi``."
            ),
            "lang": "rust",
            "load": rust_symbols,
        },
        "fortran-api": {
            "title": "Fortran",
            "intro": (
                "FORD-style comments (``!>``) and the public type-bound procedures in "
                "`dseams-fortran <https://github.com/d-SEAMS/dseams-fortran>`_. "
                "``seams_c`` is the ``bind(C)`` layer. ``seams`` is the module a program uses."
            ),
            "lang": "fortran",
            "load": fortran_symbols,
        },
        "lua-api": {
            "title": "Lua",
            "intro": (
                "Names registered on ``require(\"dseams\")`` in "
                "`yodaStruct <https://github.com/d-SEAMS/yodaStruct>`_. "
                "Each entry is a ``set_function`` name. The C column lists ``seams_*`` "
                "calls in that function when the binding reaches ``seams.h`` directly."
            ),
            "lang": "lua",
            "load": lua_symbols,
        },
        "cli-api": {
            "title": "Command line",
            "intro": (
                "Every ``seams`` subcommand and flag declared in ``src/seams_cli.cpp``. "
                ":doc:`/reference/cli` describes the text each subcommand prints. "
                "The C column lists ``seams_*`` calls in the branch that handles that subcommand."
            ),
            "lang": "bash",
            "load": cli_symbols,
        },
    }
