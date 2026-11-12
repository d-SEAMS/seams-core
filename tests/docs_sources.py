"""Each org page is in one toctree, and the book does not assume Doxygen output."""

import posixpath
import sys
from pathlib import Path

ROOT = Path("docs/orgmode")
CONF = Path("docs/source/conf.py")


def toctree_entries(text: str):
    inside = False
    for line in text.splitlines():
        if ".. toctree::" in line:
            inside = True
            continue
        if not inside:
            continue
        stripped = line.strip()
        if stripped.startswith(":"):
            continue
        if stripped == "":
            continue
        if line[:1] in (" ", "\t"):
            yield stripped
            continue
        inside = False


def docname(org_path: Path) -> str:
    rel = org_path.relative_to(ROOT).with_suffix("")
    return rel.as_posix()


def main() -> int:
    counts = {}
    for org in sorted(ROOT.rglob("*.org")):
        base = docname(org)
        parent = posixpath.dirname(base)
        for entry in toctree_entries(org.read_text(encoding="utf-8")):
            name = posixpath.normpath(posixpath.join(parent, entry))
            counts[name] = counts.get(name, 0) + 1

    problems = []
    for org in sorted(ROOT.rglob("*.org")):
        name = docname(org)
        if name == "index":
            continue
        got = counts.get(name, 0)
        if got != 1:
            problems.append(f"{name} is in {got} toctrees")

    api = counts.get("api/index", 0)
    if api != 1:
        problems.append(f"api/index is in {api} toctrees")

    conf = CONF.read_text(encoding="utf-8")
    if "if not os.path.isfile(_api_index)" not in conf:
        problems.append("conf.py does not write a C++ index when none is committed")
    if "if os.path.isfile(_xml_index)" not in conf:
        problems.append("conf.py loads doxyrest without the generated XML")

    if problems:
        sys.stderr.write("\n".join(problems) + "\n")
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
