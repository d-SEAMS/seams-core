#!/usr/bin/env python3
# Copyright (c) d-SEAMS developers
"""Positive and negative controls for check_changelog.py and check_towncrier.py.

Each case copies the changelog sources into a scratch tree, changes one
thing, and runs the checkers there as the docs suite runs them.
"""

import re
import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
FILES = (
    "CHANGELOG.md",
    "meson.build",
    "towncrier.toml",
    "docs/check_changelog.py",
    "docs/check_towncrier.py",
    "docs/towncrier.md",
    "docs/orgmode/reference/release-2.14.org",
    "scripts/project_version.py",
)
NOTES = "docs/source/reference/release-2.14-notes.rst"
PAGE = "docs/orgmode/reference/release-2.14.org"
HEADINGS = {
    "added": "Added",
    "changed": "Changed",
    "deprecated": "Deprecated",
    "removed": "Removed",
    "fixed": "Fixed",
    "security": "Security",
}
RELEASE = re.compile(r"^## \[(\d+\.\d+\.\d+)\] - \d{4}-\d{2}-\d{2}$", re.M)

sys.path.insert(0, str(ROOT / "scripts"))
from project_version import project_version  # noqa: E402


def bump(version, step):
    major, minor, patch = (int(part) for part in version.split("."))
    return "{}.{}.{}".format(major, minor, patch + step)


# The control releases follow the current one: 2.14.1 and 2.14.2 on a 2.14.0 tree.
CURRENT = project_version()
P1, P2, P3 = (bump(CURRENT, step) for step in (1, 2, 3))
FIRST_P1 = "The {} control fix is sourced from changelog.d/{}.".format(P1, P1)
SECOND_P1 = "The {} control change is sourced from changelog.d/{}.".format(P1, P1)
ONLY_P2 = "The {} control fix is sourced from changelog.d/{}.".format(P2, P2)


def section(version, date, entries):
    """CHANGELOG.md text towncrier renders for (type, slug, body) entries."""
    lines = ["## [{}] - {}".format(version, date), ""]
    for typ, title in HEADINGS.items():
        bodies = sorted(body for kind, _slug, body in entries if kind == typ)
        if not bodies:
            continue
        lines.extend(["### {}".format(title), ""])
        lines.extend("- {}".format(body) for body in bodies)
        lines.append("")
    return "\n".join(lines) + "\n"


class Tree:
    def __init__(self, path):
        self.path = path
        for name in FILES:
            target = path / name
            target.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy(ROOT / name, target)
        shutil.copytree(ROOT / "changelog.d", path / "changelog.d")

    def read(self, name):
        return (self.path / name).read_text(encoding="utf-8")

    def write(self, name, text):
        target = self.path / name
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_text(text, encoding="utf-8")

    def run(self, script, *args):
        done = subprocess.run(
            [sys.executable, str(self.path / "docs" / script), *args],
            cwd=str(self.path),
            check=False,
            capture_output=True,
            text=True,
        )
        return done.returncode, done.stdout + done.stderr

    def span(self, version):
        text = self.read("CHANGELOG.md")
        matches = list(RELEASE.finditer(text))
        for index, match in enumerate(matches):
            if match.group(1) == version:
                end = matches[index + 1].start() if index + 1 < len(matches) else len(text)
                return text, match.start(), end
        raise KeyError(version)

    def edit_section(self, version, change):
        text, start, end = self.span(version)
        self.write("CHANGELOG.md", text[:start] + change(text[start:end]) + text[end:])

    def set_version(self, version):
        meson = re.sub(
            r"^(\s*version:\s*')[^']+'",
            r"\g<1>{}'".format(version),
            self.read("meson.build"),
            count=1,
            flags=re.M,
        )
        self.write("meson.build", meson)
        page = re.sub(
            r"^#\+TITLE: .*$",
            "#+TITLE: {} release notes".format(version),
            self.read(PAGE),
            count=1,
            flags=re.M,
        )
        self.write(PAGE, page)

    def add_release(self, version, date, entries, changelog=True):
        for typ, slug, body in entries:
            self.write("changelog.d/{}/+{}.{}.md".format(version, slug, typ), body + "\n")
        if changelog:
            text = self.read("CHANGELOG.md")
            at = RELEASE.search(text).start()
            self.write("CHANGELOG.md", text[:at] + section(version, date, entries) + "\n" + text[at:])
        self.set_version(version)

    def patch_releases(self):
        self.add_release(
            P1,
            "2026-12-01",
            [("fixed", "control-fix", FIRST_P1), ("changed", "control-change", SECOND_P1)],
        )
        self.add_release(P2, "2026-12-02", [("fixed", "control-fix", ONLY_P2)])


def single_lines(text, heading=None):
    """Indexes of one-line bullets, optionally under one ### heading."""
    lines = text.splitlines(keepends=True)
    current = None
    found = []
    for index, line in enumerate(lines):
        if line.startswith("### "):
            current = line[4:].strip()
            continue
        nxt = lines[index + 1] if index + 1 < len(lines) else ""
        if line.startswith("- ") and not nxt.startswith("  ") and heading in (None, current):
            found.append(index)
    return lines, found


def drop_first(text):
    lines, found = single_lines(text)
    del lines[found[0]]
    return "".join(lines)


def edit_first(text):
    lines, found = single_lines(text)
    lines[found[0]] = lines[found[0]].rstrip("\n") + " stale\n"
    return "".join(lines)


def swap_first_two(text):
    lines, found = single_lines(text, "Added")
    a, b = found[0], found[1]
    assert b == a + 1
    lines[a], lines[b] = lines[b], lines[a]
    return "".join(lines)


def append_entry(text):
    body = text.rstrip("\n")
    return body + "\n- Extra tracked entry that is not in the fragments.\n\n"


def move_under(text, source, destination):
    lines, found = single_lines(text, source)
    entry = lines.pop(found[0])
    at = next(i for i, line in enumerate(lines) if line == "### {}\n".format(destination))
    lines.insert(at + 2, entry)
    return "".join(lines)


class ChangelogChecks(unittest.TestCase):
    def setUp(self):
        self.scratch = tempfile.TemporaryDirectory(prefix="seams-changelog-controls-")
        self.tree = Tree(Path(self.scratch.name))

    def tearDown(self):
        self.scratch.cleanup()

    def passes(self, script, *args):
        code, output = self.tree.run(script, *args)
        self.assertEqual(code, 0, output)
        return output

    def fails(self, script, *needles):
        code, output = self.tree.run(script, *(["--check"] if script == "check_changelog.py" else []))
        self.assertEqual(code, 1, output)
        for needle in needles:
            self.assertIn(needle, output)
        return output

    # Clean trees.

    def test_clean_tree_passes(self):
        self.passes("check_changelog.py", "--check")
        self.passes("check_towncrier.py")

    def test_check_writes_nothing_into_the_tree(self):
        self.passes("check_changelog.py", "--check")
        self.assertFalse((self.tree.path / NOTES).exists())

    def test_unknown_arguments_are_usage_errors(self):
        for args in (["--output", "notes.rst"], ["--selftest"], ["--check", "--check"]):
            with self.subTest(args=args):
                code, output = self.tree.run("check_changelog.py", *args)
                self.assertEqual(code, 2, output)
                self.assertIn("usage:", output)

    def test_valid_patch_releases_pass(self):
        self.tree.patch_releases()
        self.assertIn("{} 1, {} 2, ".format(P2, P1), self.passes("check_towncrier.py"))
        self.passes("check_changelog.py", "--check")

    # The 2.14.0 section built from the top-level fragments.

    def test_missing_entry(self):
        self.tree.edit_section("2.14.0", drop_first)
        self.fails("check_towncrier.py", "2.14.0 is missing an entry under Added")

    def test_extra_entry_after_the_draft(self):
        _text, start, end = self.tree.span("2.14.0")
        draft = _text[start:end].strip()
        self.tree.edit_section("2.14.0", append_entry)
        self.assertIn(draft, self.tree.read("CHANGELOG.md"))
        self.fails("check_towncrier.py", "2.14.0 has an extra entry under Fixed: Extra tracked entry")

    def test_edited_entry(self):
        self.tree.edit_section("2.14.0", edit_first)
        self.fails("check_towncrier.py", "2.14.0 has a stale or edited entry under Added")

    def test_entries_out_of_towncrier_order(self):
        self.tree.edit_section("2.14.0", swap_first_two)
        self.fails("check_towncrier.py", "2.14.0 entries are not in towncrier order")

    def test_entry_under_the_wrong_type(self):
        self.tree.edit_section("2.14.0", lambda text: move_under(text, "Fixed", "Added"))
        self.fails(
            "check_towncrier.py",
            "entry is in the wrong section (found in 2.14.0 Added, expected in 2.14.0 Fixed)",
        )

    def test_entry_in_an_older_release(self):
        text, start, end = self.tree.span("2.14.0")
        lines, found = single_lines(text[start:end], "Added")
        entry = lines.pop(found[0])
        older = text[end:]
        heading = re.search(r"^### (\w+)\n\n", older, re.M)
        older = older[: heading.end()] + entry + older[heading.end() :]
        self.tree.write("CHANGELOG.md", text[:start] + "".join(lines) + older)
        self.fails(
            "check_towncrier.py",
            "entry is in the wrong section (found in 2.13.0 {}, expected in 2.14.0 Added)".format(
                heading.group(1)
            ),
        )

    def test_removed_fragment(self):
        (self.tree.path / "changelog.d" / "+mpi-measured-cut.fixed.md").unlink()
        self.fails("check_towncrier.py", "2.14.0 has an extra entry under Fixed: An MPI build")
        self.passes("check_changelog.py", "--check")

    def test_removed_final_fragment(self):
        text, start, end = self.tree.span("2.14.0")
        last = [line for line in text[start:end].splitlines() if line.startswith("- ")][-1][2:]
        owner = [
            path
            for path in (self.tree.path / "changelog.d").glob("*.md")
            if path.read_text(encoding="utf-8").strip() == last
        ]
        self.assertEqual(len(owner), 1)
        owner[0].unlink()
        self.fails("check_towncrier.py", "2.14.0 has an extra entry under Fixed: " + last[:40])

    def test_edited_fragment(self):
        path = self.tree.path / "changelog.d" / "+mpi-measured-cut.fixed.md"
        path.write_text(path.read_text(encoding="utf-8").strip() + " stale\n", encoding="utf-8")
        self.fails("check_towncrier.py", "2.14.0 has a stale or edited entry under Fixed: An MPI build")

    def test_duplicate_fragment(self):
        body = self.tree.read("changelog.d/+mpi-measured-cut.fixed.md")
        self.tree.write("changelog.d/+duplicate-control.fixed.md", body)
        self.fails("check_changelog.py", "repeats")

    def test_too_few_fragments(self):
        for path in sorted((self.tree.path / "changelog.d").glob("+*.md"))[:200]:
            path.unlink()
        self.fails("check_changelog.py", "expected at least 100 towncrier fragments")
        self.fails("check_towncrier.py", "changelog.d lost fragments")

    # Patch releases in changelog.d/<X.Y.Z>.

    def test_missing_entry_of_the_first_patch_release(self):
        self.tree.patch_releases()
        self.tree.write("changelog.d/{}/+control-extra.fixed.md".format(P1), "A fix nobody listed.\n")
        self.fails("check_towncrier.py", "{} is missing an entry under Fixed: A fix nobody".format(P1))

    def test_missing_entry_of_the_second_patch_release(self):
        self.tree.patch_releases()
        self.tree.write("changelog.d/{}/+control-extra.fixed.md".format(P2), "A fix nobody listed.\n")
        self.fails("check_towncrier.py", "{} is missing an entry under Fixed: A fix nobody".format(P2))

    def test_patch_entry_under_the_wrong_type(self):
        self.tree.patch_releases()
        self.tree.edit_section(P2, lambda text: text.replace("### Fixed", "### Added"))
        self.fails(
            "check_towncrier.py",
            "entry is in the wrong section (found in {0} Added, expected in {0} Fixed)".format(P2),
        )

    def test_patch_entry_in_the_wrong_release(self):
        self.tree.patch_releases()
        first = "- {}\n".format(FIRST_P1)
        only = "- {}\n".format(ONLY_P2)
        self.tree.edit_section(P1, lambda text: text.replace(first, ""))
        self.tree.edit_section(P2, lambda text: text.replace(only, first + only))
        self.fails(
            "check_towncrier.py",
            "entry is in the wrong section (found in {} Fixed, expected in {} Fixed)".format(P2, P1),
        )
        self.fails(
            "check_changelog.py",
            "CHANGELOG.md [{}] is missing the towncrier bullet for: +control-fix.fixed.md".format(P1),
        )

    def test_patch_entry_with_stale_continuation(self):
        self.tree.patch_releases()
        first = "- {}\n".format(FIRST_P1)
        self.tree.edit_section(
            P1,
            lambda text: text.replace(first, first + "  stale continuation\n", 1),
        )
        self.fails(
            "check_changelog.py",
            "CHANGELOG.md [{}] is missing the towncrier bullet for: +control-fix.fixed.md".format(P1),
        )

    def test_patch_section_without_fragments(self):
        self.tree.patch_releases()
        shutil.rmtree(self.tree.path / "changelog.d" / P2)
        self.fails("check_towncrier.py", "CHANGELOG.md [{0}] has no changelog.d/{0} fragments".format(P2))

    def test_fragments_without_a_section(self):
        self.tree.patch_releases()
        self.tree.add_release(
            P3, "2026-12-03", [("fixed", "control-fix", "An unreleased fix.")], changelog=False
        )
        self.fails("check_towncrier.py", "CHANGELOG.md has no [{}] section for its fragments".format(P3))

    def test_repeated_heading(self):
        self.tree.patch_releases()
        _text, start, end = self.tree.span(P2)
        self.tree.edit_section(P1, lambda text: _text[start:end] + text)
        self.fails("check_towncrier.py", "CHANGELOG.md repeats the [{}] heading".format(P2))

    def test_headings_not_newest_first(self):
        self.tree.patch_releases()
        text, start, end = self.tree.span(P2)
        moved = text[:start] + text[end:]
        at = moved.index("## [2.14.0]")
        self.tree.write("CHANGELOG.md", moved[:at] + text[start:end] + moved[at:])
        self.fails("check_towncrier.py", "CHANGELOG.md releases are not newest first")

    def test_newest_heading_is_not_the_project_version(self):
        self.tree.patch_releases()
        self.tree.set_version(P1)
        self.fails(
            "check_changelog.py",
            "newest release heading is {}, meson.build version is {}".format(P2, P1),
        )

    def test_page_title_is_not_the_project_version(self):
        self.tree.patch_releases()
        page = self.tree.read(PAGE).replace("#+TITLE: " + P2, "#+TITLE: " + P1)
        self.tree.write(PAGE, page)
        self.fails(
            "check_changelog.py",
            "title is '{} release notes', expected '{} release notes'".format(P1, P2),
        )

    def test_bad_name_in_a_patch_release(self):
        self.tree.patch_releases()
        self.tree.write("changelog.d/{}/control.fixed.md".format(P2), "A fix without the prefix.\n")
        self.fails("check_changelog.py", "{}/control.fixed.md does not match fragment naming".format(P2))

    def test_empty_fragment_in_a_patch_release(self):
        self.tree.patch_releases()
        self.tree.write("changelog.d/{}/+control-empty.fixed.md".format(P2), "\n")
        self.fails("check_changelog.py", "{}/+control-empty.fixed.md is empty".format(P2))

    def test_repeat_across_releases(self):
        self.tree.patch_releases()
        body = self.tree.read("changelog.d/+mpi-measured-cut.fixed.md")
        self.tree.write("changelog.d/{}/+control-repeat.fixed.md".format(P2), body)
        self.fails(
            "check_changelog.py",
            "{}/+control-repeat.fixed.md repeats +mpi-measured-cut.fixed.md".format(P2),
        )

    def test_empty_patch_release_directory(self):
        (self.tree.path / "changelog.d" / P1).mkdir()
        self.fails("check_changelog.py", "changelog.d/{} has no fragments".format(P1))

    def test_directory_outside_the_series(self):
        self.tree.write("changelog.d/2.13.5/+control.fixed.md", "A fix for another series.\n")
        self.fails("check_changelog.py", "changelog.d/2.13.5 is not a 2.14 patch release")
        self.fails("check_towncrier.py", "CHANGELOG.md has no [2.13.5] section for its fragments")

    def test_directory_that_is_not_a_release(self):
        self.tree.write("changelog.d/notes/+control.fixed.md", "A fix in a stray directory.\n")
        self.fails("check_changelog.py", "changelog.d/notes is not a release directory")
        self.fails("check_towncrier.py", "changelog.d/notes is not a release directory")

    # The rendering check_changelog.py validates.

    def test_rendering_that_does_not_match_the_fragments(self):
        self.tree.patch_releases()
        probe = (
            "import sys\n"
            "sys.path.insert(0, 'docs')\n"
            "import check_changelog as c\n"
            "assert c.main(['--check']) == 0\n"
            "errors, fragments = c.load()\n"
            "good = c.render(fragments)\n"
            "version = c.project_version()\n"
            "assert c.check(fragments, good, version) == []\n"
            "first = good.split('\\n- ', 1)[1].split('\\n', 1)[0]\n"
            "cases = {\n"
            "    'edited': good.replace(first, first + ' mismatch', 1),\n"
            "    'dropped': good.replace('- ' + first + '\\n', '', 1),\n"
            "    'repeated': good.replace('- ' + first + '\\n', '- ' + first + '\\n\\n- ' + first + '\\n', 1),\n"
            "    'unsourced': good.replace('{0}\\n------\\n', '{1}\\n------\\n', 1),\n"
            "}\n"
            "for name, text in cases.items():\n"
            "    assert text != good, name\n"
            "    print(name, '|'.join(c.check(fragments, text, version)))\n"
        )
        done = subprocess.run(
            [sys.executable, "-c", probe.replace("{0}", P2).replace("{1}", P3)],
            cwd=str(self.tree.path),
            check=False,
            capture_output=True,
            text=True,
        )
        self.assertEqual(done.returncode, 0, done.stdout + done.stderr)
        results = dict(line.split(" ", 1) for line in done.stdout.splitlines() if " " in line)
        self.assertIn("orphan fragment", results["edited"])
        self.assertIn("note is not a fragment", results["edited"])
        self.assertIn("orphan fragment", results["dropped"])
        self.assertIn("duplicate note", results["repeated"])
        self.assertIn("release headings are", results["unsourced"])


if __name__ == "__main__":
    unittest.main()
