"""seams under mpiexec -n 2 prints what the serial run prints.

usage: mpi_frames.py MPIEXEC SEAMS DUMP CON
"""

import os
import signal
import subprocess
import sys
import tempfile

mpiexec, seams, dump, con = sys.argv[1:5]
# Colour allowed but not forced, so a rank writing it into piped output
# differs from serial
base = dict(os.environ, TERM="xterm-256color")
for name in ("NO_COLOR", "FORCE_COLOR", "CLICOLOR", "CLICOLOR_FORCE",
             "SEAMS_MPI"):
    base.pop(name, None)
# A bare run that started MPI anyway would fail on this transport
bare = dict(base, OMPI_MCA_pml="no_such_pml")
failed = []


def run(cmd, env):
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                         text=True, env=env, start_new_session=True)
    try:
        out, err = p.communicate(timeout=90)
    except subprocess.TimeoutExpired:
        # A hung rank must not outlive the test: mpiexec passes SIGTERM on
        # to its ranks, which its launcher daemon keeps out of this group
        p.terminate()
        try:
            p.communicate(timeout=20)
        except subprocess.TimeoutExpired:
            os.killpg(p.pid, signal.SIGKILL)
            p.communicate()
        raise
    return subprocess.CompletedProcess(cmd, p.returncode, out, err)


def check(name, ok):
    if not ok:
        failed.append(name)


def same(args, records, rc=0, extra=None):
    serial = run([seams, "--format", "json"] + args, dict(bare, **(extra or {})))
    mpi = run([mpiexec, "-n", "2", seams, "--format", "json"] + args,
              dict(base, **(extra or {})))
    check(" ".join(args), serial.returncode == rc and
          (mpi.returncode == 0) == (rc == 0) and serial.stdout == mpi.stdout and
          serial.stdout.count('{"schema"') == records)


# A Steinhardt atom split still spanning both ranks pairs different frames
# and hangs
same(["--frame", "1", "--last", "11", "steinhardt", dump], 11)
same(["--frame", "1", "--last", "20", "cages", dump], 11)
# read prints the frame it loaded, so a rank reading the wrong frame shows;
# frame 0 starts at the first
same(["--frame", "0", "--last", "20", "read", dump], 11)
# The CON reader stops at its first missing frame and stays on rank 0
same(["--frame", "1", "--last", "8", "read", con], 2)
same(["--frame", "3", "steinhardt", dump], 1)
# Every frame fails, and each status comes back through the gather
same(["--strict-input", "--type", "99", "--frame", "1", "--last", "11", "read",
      dump], 11, rc=2)
same(["--frame", "1", "--last", "4", "read", dump], 4, extra={"FORCE_COLOR": "1"})

with tempfile.TemporaryDirectory() as tmp:
    out = os.path.join(tmp, "cages.dump")
    p = run([mpiexec, "-n", "2", seams, "--frame", "1", "--last", "4",
             "--per-atom", out, "cages", dump], base)
    check("--per-atom on two ranks", p.returncode != 0 and not os.path.exists(out))

# Two commands, or two ranges, in one launch are refused
for name, a, b in (
        ("two commands", [seams, "--frame", "1", "--last", "4", "read", dump],
         [seams, "--frame", "1", "--last", "4", "steinhardt", dump]),
        ("two ranges", ["env", "SEAMS_LAST=11", seams, "--frame", "1", "read", dump],
         ["env", "SEAMS_LAST=5", seams, "--frame", "1", "read", dump])):
    p = run([mpiexec, "-n", "1"] + a + [":", "-n", "1"] + b, base)
    check("refuse " + name, p.returncode == 2 and "nop" not in p.stdout)

# SEAMS_MPI=0 keeps MPI off: each process prints the whole range, on the
# terminal Open MPI gives it
args = ["--format", "json", "--frame", "1", "--last", "4", "read", dump]
serial = run([seams] + args, dict(bare, NO_COLOR="1"))
p = run([mpiexec, "-n", "2", "env", "SEAMS_MPI=0", seams] + args,
        dict(base, NO_COLOR="1"))
check("SEAMS_MPI=0", p.returncode == 0 and
      sorted(p.stdout.splitlines()) == sorted(serial.stdout.splitlines() * 2))

if failed:
    print("failed:", *failed, sep="\n  ")
sys.exit(1 if failed else 0)
