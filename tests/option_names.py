"""CLI flags and Meson project options do not share or shadow a name."""

import re
import subprocess
import sys
from pathlib import Path


def long_names(help_text: str):
    names = []
    for line in help_text.splitlines():
        stripped = line.strip()
        if not stripped.startswith("--"):
            continue
        token = stripped.split()[0].rstrip(",")
        if token.startswith("--"):
            names.append(token)
    return names


def declared_options(path: Path):
    if not path.is_file():
        return set()
    return set(re.findall(r"option\(\s*'([^']+)'", path.read_text(encoding="utf-8")))


def configured_project_options(build_dir: str):
    proc = subprocess.run(
        ["meson", "configure", build_dir],
        check=False,
        text=True,
        capture_output=True,
    )
    if proc.returncode != 0:
        sys.stderr.write(proc.stdout)
        sys.stderr.write(proc.stderr)
        raise SystemExit(proc.returncode)
    current = None
    found = {}
    in_project = False
    name_re = re.compile(r"^  ([A-Za-z_][A-Za-z0-9_]*)\s{2,}")
    for line in proc.stdout.splitlines():
        if line.startswith("Subproject "):
            current = line.split("Subproject ", 1)[1].rstrip(":").strip()
            found.setdefault(current, set())
            in_project = False
            continue
        if current is None:
            continue
        if line.startswith("  Project options"):
            in_project = True
            continue
        if not in_project:
            continue
        if not line.strip() or line.strip().startswith("-"):
            continue
        match = name_re.match(line)
        if match:
            found[current].add(match.group(1))
            continue
        in_project = False
    return found


def argum_ambiguous(names, token):
    """True when Argum cannot bind token to one long option.

    findMatchOrMatchingPrefixRange returns an exact key by itself, then
    stops. A longer option that starts with that key is not a second
    candidate. A token that is not itself an option is ambiguous when
    two or more options start with it.
    """
    if token in names:
        return False
    return sum(name.startswith(token) for name in names) > 1


def main() -> int:
    if len(sys.argv) != 3:
        sys.stderr.write("usage: option_names.py SEAMS BUILDDIR\n")
        return 2
    seams, build_dir = sys.argv[1], sys.argv[2]
    proc = subprocess.run([seams, "--help"], check=False, text=True, capture_output=True)
    if proc.returncode != 0:
        sys.stderr.write(proc.stderr)
        return proc.returncode
    names = long_names(proc.stdout)
    if len(names) != len(set(names)):
        sys.stderr.write("duplicate long options in --help\n")
        return 1
    # Exact --type and --types both bind. Only a non-option prefix such
    # as --ty is ambiguous, and that token is not a registered name.
    sample = {"--type", "--types", "--cutoff"}
    if (argum_ambiguous(sample, "--type") or argum_ambiguous(sample, "--types")
            or not argum_ambiguous(sample, "--ty")
            or argum_ambiguous(sample, "--cu")):
        sys.stderr.write("Argum prefix rule mismatch\n")
        return 1
    known = set(names)
    for left in names:
        for right in names:
            if left != right and right.startswith(left) and argum_ambiguous(known, left):
                sys.stderr.write(f"{left} is ambiguous with {right}\n")
                return 1

    for name in ("gpulite", "vesin"):
        patch = Path("subprojects/packagefiles") / name / "meson_options.txt"
        if not patch.is_file():
            sys.stderr.write(f"missing {patch}\n")
            return 1
    listed = configured_project_options(build_dir)
    root = Path("subprojects")
    # nanobind and its header map are optional wraps. They are not seams
    # projects and do not carry this tree's option file.
    third_party = {"nanobind", "robin-map"}
    for project, options in sorted(listed.items()):
        if project in third_party:
            continue
        src = root / project / "meson_options.txt"
        if not src.is_file():
            sys.stderr.write(f"{project} has no meson_options.txt\n")
            return 1
        extra = options - declared_options(src)
        if extra:
            sys.stderr.write(
                f"{project} lists options it does not declare: {sorted(extra)}\n"
            )
            return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
