#!/usr/bin/env python3
"""Each engine page names the flags the source actually accepts.

The pages are docs/source/engines. The flags are read from the engine
sources, not from a second copy of the documentation. Pass a build
directory to also run the minimal example for every engine.
"""

import ast
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
PAGES = ROOT / "docs" / "source" / "engines"


def _text(path):
    return (ROOT / path).read_text(encoding="utf-8")


def gromacs_flags():
    return set(re.findall(r'DoubleOption\("([^"]+)"\)', _text("engines/gromacs/seams_ice.cpp")))


def lammps_flags():
    source = _text("engines/lammps/compute_dseams.cpp")
    flags = set(re.findall(r'strcmp\(arg\[iarg\], "([^"]+)"\)', source))
    if "pf/callback" not in _text("engines/lammps/in.dseams"):
        raise SystemExit("engines/lammps/in.dseams does not install fix external pf/callback")
    flags.add("pf/callback")
    return flags


def plumed_flags():
    source = _text("engines/plumed/SeamsColvar.cpp")
    flags = set(re.findall(r'keys\.add\("[^"]+", "([^"]+)"', source))
    flags.update(re.findall(r'addOutputComponent\("([^"]+)"', source))
    return flags


def openmm_flags():
    source = _text("adapters/python/dseams_insitu/openmm.py")
    match = re.search(r"def __init__\(self, (.*?)\):", source, re.S)
    if match is None:
        raise SystemExit("DseamsReporter has no constructor")
    names = []
    for part in match.group(1).split(","):
        name = part.strip().split("=")[0].strip()
        if name and name != "self":
            names.append(name)
    if "DseamsReporter" not in source:
        raise SystemExit("openmm adapter is missing DseamsReporter")
    return set(names)


def insitu_flags():
    source = _text("adapters/python/dseams_insitu/engine.py")
    module = ast.parse(source)
    found = {}
    for node in module.body:
        if isinstance(node, ast.ClassDef) and node.name == "InSitu":
            for item in node.body:
                if isinstance(item, ast.FunctionDef) and item.name in {"__init__", "observe"}:
                    found[item.name] = [arg.arg for arg in item.args.args if arg.arg != "self"]
    if set(found) != {"__init__", "observe"}:
        raise SystemExit("InSitu is missing __init__ or observe")
    return set(found["__init__"]) | set(found["observe"])


ENGINES = {
    "gromacs": gromacs_flags,
    "lammps": lammps_flags,
    "openmm": openmm_flags,
    "plumed": plumed_flags,
    "insitu": insitu_flags,
}


def page_text(page):
    """Page prose plus every file the page literal-includes."""
    text = page.read_text(encoding="utf-8")
    included = []
    for match in re.finditer(r"\{literalinclude\}\s+(\S+)", text):
        target = (page.parent / match.group(1)).resolve()
        if target.is_file():
            included.append(target.read_text(encoding="utf-8"))
    return text + "\n" + "\n".join(included)


def main(argv):
    build = Path(argv[1]) if len(argv) > 1 else None
    for name, reader in ENGINES.items():
        page = PAGES / (name + ".md")
        if not page.is_file():
            raise SystemExit(f"missing engine page {page}")
        text = page_text(page)
        flags = reader()
        if not flags:
            raise SystemExit(f"{name} source advertised no flags")
        missing = sorted(flag for flag in flags if flag not in text)
        if missing:
            raise SystemExit(f"{name} page is missing flags from the source: {', '.join(missing)}")
        print(f"{name} page names {', '.join(sorted(flags))}")
    if build is None:
        print("engine pages match the source; pass a build directory to run the examples")
        return 0
    sys.path.insert(0, str(ROOT / "docs" / "howto"))
    import test_engines

    for name in ENGINES:
        test_engines.main_engine = None
        # test_engines.main reads argv. Call the runners directly.
        {
            "gromacs": test_engines.run_gromacs,
            "lammps": test_engines.run_lammps,
            "openmm": test_engines.run_openmm,
            "plumed": test_engines.run_plumed,
            "insitu": test_engines.run_insitu,
        }[name](ROOT, build)
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
