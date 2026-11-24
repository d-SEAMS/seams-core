# Copyright (c) d-SEAMS developers
"""Configure a scratch project that takes seams-core as a Meson subproject.

The scratch project's subprojects/ links to the wraps the build under test
already unpacked (the top project's subproject directory), and the
configure does not download, so nothing is written to either source tree.
WORK_DIR is removed first and afterwards. chemfiles is off in the scratch
build: no mode reads it, and a vendored copy is compiled at configure time
on every core, starving the tests that run beside this one.

MODE with-tests-false: -Dseams-core:with_tests=false must register no test,
because a consumer runs every test its subprojects register.

MODE source-tree-hidden: -Dseams-core:with_tests=true must not register the
four source_tree_untracked_* tests. Those belong to a top-level build only.

MODE build-dirs: with -Dseams-core:with_tests=true, the scripts that take
seams-core's build directory (option_names.py, check_no_rebuild.py and the
examples/c/check_installed.py helpers) must find the parent build directory
above it, since meson configure, introspect and install accept only that.

Usage: subproject_consumer.py MODE MESON SOURCE_ROOT TOP_BUILD_ROOT SEAMS WORK_DIR

SEAMS is the seams executable for build-dirs, and is not read otherwise.
"""

import importlib.util
import json
import shutil
import subprocess
import sys
from pathlib import Path

sys.dont_write_bytecode = True
_HERE = str(Path(__file__).resolve().parent)
if _HERE not in sys.path:
    sys.path.insert(0, _HERE)
from host_tool_env import without_asan_preload  # noqa: E402

PROJECT = "seams-core"
MESON_BUILD = """project('seams-core-consumer', 'c', 'cpp')
subproject('seams-core')
"""


def top_subprojects(top_build: Path) -> Path:
    info = json.loads((top_build / "meson-info" / "meson-info.json").read_text(encoding="utf-8"))
    project = json.loads((top_build / "meson-info" / "intro-projectinfo.json").read_text(encoding="utf-8"))
    return Path(info["directories"]["source"]) / project["subproject_dir"]


def scratch_project(work: Path, source: Path, top_build: Path) -> Path:
    parent = work / "parent"
    subprojects = parent / "subprojects"
    subprojects.mkdir(parents=True)
    (parent / "meson.build").write_text(MESON_BUILD, encoding="utf-8")
    for entry in sorted(top_subprojects(top_build).iterdir()):
        if entry.name in (PROJECT, "packagecache"):
            continue
        if entry.suffix == ".wrap":
            shutil.copy(entry, subprojects / entry.name)
        elif entry.is_dir():
            (subprojects / entry.name).symlink_to(entry.resolve(), target_is_directory=True)
    (subprojects / PROJECT).symlink_to(source.resolve(), target_is_directory=True)
    return parent


def run(cmd, cwd=None) -> subprocess.CompletedProcess:
    return subprocess.run(
        [str(c) for c in cmd],
        cwd=cwd,
        capture_output=True,
        text=True,
        env=without_asan_preload(),
    )


def helpers(source: Path):
    spec = importlib.util.spec_from_file_location("quickstart_installed", source / "examples" / "c" / "check_installed.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


_SOURCE_TREE_TESTS = (
    "source_tree_untracked_snapshot",
    "source_tree_untracked_check",
    "source_tree_untracked_selftest",
    "source_tree_untracked_parallel",
)


def listed_tests(meson: str, build: Path) -> tuple[list[str], list[dict]]:
    listed = run([meson, "introspect", "--tests", build])
    if listed.returncode != 0:
        return [listed.stdout + listed.stderr + "meson introspect --tests failed"], []
    return [], json.loads(listed.stdout)


def with_tests_false(meson: str, build: Path) -> list[str]:
    errors, tests = listed_tests(meson, build)
    if errors:
        return errors
    return [f"registered with -Dseams-core:with_tests=false: {t['name']}" for t in tests]


def source_tree_hidden(meson: str, build: Path) -> list[str]:
    errors, tests = listed_tests(meson, build)
    if errors:
        return errors
    names = [t["name"] for t in tests]
    hidden = [name for name in names if any(name.endswith(stem) for stem in _SOURCE_TREE_TESTS)]
    if hidden:
        return [f"registered in a subproject: {name}" for name in hidden]
    if not names:
        return ["-Dseams-core:with_tests=true registered no tests"]
    return []


def build_dirs(meson: str, source: Path, seams: str, build: Path) -> list[str]:
    errors = []
    project_build = build / "subprojects" / PROJECT
    python = sys.executable
    for cmd in (
        [python, source / "tests" / "option_names.py", seams, project_build],
        [python, source / "tests" / "check_no_rebuild.py", source, project_build],
    ):
        done = run(cmd, cwd=source)
        if done.returncode != 0:
            errors.append(f"{Path(cmd[1]).name} exited {done.returncode}:\n{done.stdout}{done.stderr}")
    help_mod = helpers(source)
    want = json.loads(run([meson, "introspect", "--projectinfo", build]).stdout)
    version = next(s["version"] for s in want["subprojects"] if s["name"] == PROJECT)
    try:
        root = help_mod.meson_root(project_build)
        found = help_mod.project_version(project_build)
    except (AttributeError, subprocess.CalledProcessError, SystemExit) as exc:
        errors.append(f"examples/c/check_installed.py on {project_build}: {exc!r}")
        return errors
    if root != build:
        errors.append(f"meson_root({project_build}) is {root}, not {build}")
    if found != version:
        errors.append(f"project_version gave {found}, not {version}")
    if not (root / "meson-uninstalled" / "seams-uninstalled.pc").is_file():
        errors.append("meson-uninstalled/seams-uninstalled.pc is not in the parent build directory")
    return errors


def main(argv: list[str]) -> int:
    if len(argv) != 7 or argv[1] not in ("with-tests-false", "source-tree-hidden", "build-dirs"):
        print(__doc__, file=sys.stderr)
        return 2
    mode, meson = argv[1], argv[2]
    source, top_build, seams, work = Path(argv[3]), Path(argv[4]), argv[5], Path(argv[6])
    with_tests = "false" if mode == "with-tests-false" else "true"
    shutil.rmtree(work, ignore_errors=True)
    try:
        parent = scratch_project(work, source, top_build)
        build = work / "build"
        setup = run([
            meson, "setup", build, parent, f"-D{PROJECT}:with_tests={with_tests}",
            f"-D{PROJECT}:with_chemfiles=disabled", "--wrap-mode=nodownload",
        ])
        if setup.returncode != 0:
            print(setup.stdout + setup.stderr, file=sys.stderr)
            print(f"meson setup -D{PROJECT}:with_tests={with_tests} exited {setup.returncode}", file=sys.stderr)
            return 1
        if mode == "with-tests-false":
            errors = with_tests_false(meson, build)
        elif mode == "source-tree-hidden":
            errors = source_tree_hidden(meson, build)
        else:
            errors = build_dirs(meson, source, seams, build.resolve())
    finally:
        shutil.rmtree(work, ignore_errors=True)
    for error in errors:
        print(error, file=sys.stderr)
    if errors:
        return 1
    if mode == "with-tests-false":
        print(f"-D{PROJECT}:with_tests=false registers no tests in a consumer")
    elif mode == "source-tree-hidden":
        print("a subproject with tests on does not register the source-tree checks")
    else:
        print("build-directory scripts find the parent build directory of a subproject")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
