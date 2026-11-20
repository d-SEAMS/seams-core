"""`seams read` fails on a missing or unreadable structure file.

usage: read_unreadable.py SEAMS EXT

EXT is pdb, gro, dcd, or con. For that one format the script runs
`seams read` on three paths:

  * a path that does not exist
  * an empty file
  * a file whose mode is 000

Each run must exit 2, write nothing to stdout, and print this stderr line:

  cannot read <path>

A file the reader can open and that contains atoms is not part of this
test. The check does not skip.
"""

import os
import stat
import subprocess
import sys
import tempfile


EXTENSIONS = ("pdb", "gro", "dcd", "con")


def check(seams, path):
    proc = subprocess.run([seams, "read", path], capture_output=True, text=True)
    expected = "cannot read " + path
    lines = proc.stderr.splitlines()
    if proc.returncode == 2 and proc.stdout == "" and expected in lines:
        print("ok " + os.path.basename(path) + " exit " + str(proc.returncode))
        return None
    return (
        "path " + path + "\n"
        "exit " + str(proc.returncode) + "\n"
        "stdout " + repr(proc.stdout) + "\n"
        "stderr " + repr(proc.stderr) + "\n"
        "expected exit 2, empty stdout, and a stderr line " + repr(expected)
    )


def main():
    if len(sys.argv) != 3:
        print("usage: read_unreadable.py SEAMS EXT", file=sys.stderr)
        return 2
    seams = sys.argv[1]
    ext = sys.argv[2]
    if ext not in EXTENSIONS:
        print("extension " + repr(ext) + " is not one of " + str(EXTENSIONS), file=sys.stderr)
        return 2
    if not os.path.isfile(seams):
        print("seams binary is missing: " + seams, file=sys.stderr)
        return 2

    failures = []
    tmp = tempfile.mkdtemp(prefix="read-unreadable-" + ext + "-")
    try:
        missing = os.path.join(tmp, "missing." + ext)
        empty = os.path.join(tmp, "empty." + ext)
        with open(empty, "wb") as handle:
            handle.write(b"")
        blocked = os.path.join(tmp, "blocked." + ext)
        with open(blocked, "wb") as handle:
            handle.write(b"blocked")
        os.chmod(blocked, 0)
        try:
            for path in (missing, empty, blocked):
                detail = check(seams, path)
                if detail is not None:
                    print(detail, file=sys.stderr)
                    failures.append(os.path.basename(path))
        finally:
            os.chmod(blocked, stat.S_IRUSR | stat.S_IWUSR)
    finally:
        for name in os.listdir(tmp):
            os.remove(os.path.join(tmp, name))
        os.rmdir(tmp)
    if failures:
        print("failed: " + ", ".join(failures), file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
