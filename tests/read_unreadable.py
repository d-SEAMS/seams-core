"""`seams read` fails on a missing or unreadable structure file.

usage: read_unreadable.py SEAMS EXT

EXT is pdb, gro, dcd, or con. For that one format the script runs
`seams read` and `seams --frame 1 --last 2 read` on four paths:

  * a path that does not exist
  * an empty file
  * a file whose mode is 000, which this process cannot read
  * a directory

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


def check(seams, args, path):
    command = [seams, *args, path]
    try:
        proc = subprocess.run(command, capture_output=True, text=True, timeout=30)
    except subprocess.TimeoutExpired as exc:
        return (
            "command " + " ".join(command) + "\n"
            "timed out\n"
            "stdout " + repr(exc.stdout) + "\n"
            "stderr " + repr(exc.stderr)
        )
    expected = "cannot read " + path
    lines = proc.stderr.splitlines()
    shown = " ".join(args + [os.path.basename(path)])
    if proc.returncode == 2 and proc.stdout == "" and expected in lines:
        print("ok " + shown + " exit " + str(proc.returncode))
        return None
    return (
        "command " + " ".join(command) + "\n"
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
        directory = os.path.join(tmp, "directory." + ext)
        os.mkdir(directory)
        try:
            with open(blocked, "rb") as handle:
                handle.read(1)
        except PermissionError:
            denied = True
        else:
            denied = False
        if not denied:
            print("mode 000 file is still readable: " + blocked, file=sys.stderr)
            failures.append(os.path.basename(blocked))
        try:
            commands = (
                ["read"],
                ["--frame", "1", "--last", "2", "read"],
            )
            for path in (missing, empty, blocked, directory):
                for args in commands:
                    detail = check(seams, args, path)
                    if detail is not None:
                        print(detail, file=sys.stderr)
                        failures.append(" ".join(args + [os.path.basename(path)]))
        finally:
            os.chmod(blocked, stat.S_IRUSR | stat.S_IWUSR)
    finally:
        for name in os.listdir(tmp):
            path = os.path.join(tmp, name)
            if os.path.isdir(path) and not os.path.islink(path):
                os.rmdir(path)
            else:
                os.remove(path)
        os.rmdir(tmp)
    if failures:
        print("failed: " + ", ".join(failures), file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
