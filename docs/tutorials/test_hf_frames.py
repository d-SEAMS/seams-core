#!/usr/bin/env python3
"""Nucleation, interface, and brine tutorials.

``--small`` is the CI mode. It runs the shipped frames:

* ``input/traj/interface_f0.lammpstrj``
* ``input/traj/nucleation_three.lammpstrj``
* ``input/traj/nacl_front.lammpstrj``

Without ``--small``, each tutorial tries a pinned Hugging Face revision
and uses the shipped frame when that fetch does not return the file.
Brine is ``d-SEAMS/brine-front-230K`` at
``4828e33ecf419c70643b2d466bfa4f283b4138dc``. Interface is
``d-SEAMS/interface-frames`` and nucleation is
``d-SEAMS/nucleation-frames``. Those two revisions are empty until the
datasets are uploaded, so the fetch falls back to the in-tree dumps.
``scripts/hf_dataset_prepare.py`` packs those dumps and does not upload them.

    python docs/tutorials/test_hf_frames.py --small ROOT SEAMS LIBDIR
"""

import hashlib
import sys
import urllib.error
import urllib.request
from pathlib import Path

PIN_REPO = "d-SEAMS/brine-front-230K"
PIN_REV = "4828e33ecf419c70643b2d466bfa4f283b4138dc"
PIN_FILE = "data/frames/T230_m20_r0_tend.xyz"
PIN_SHA256 = "dcbc05f81f321b45195d259cd4c22e69340f9c1ae57d93c47a74c995c58dd05d"
PIN_URL = "https://huggingface.co/datasets/{}/resolve/{}/{}".format(PIN_REPO, PIN_REV, PIN_FILE)

# Empty until the packed datasets are uploaded. A set revision is fetched;
# any failure uses the shipped dump.
INTERFACE_REPO = "d-SEAMS/interface-frames"
INTERFACE_REV = ""
INTERFACE_REMOTE = "data/frames/interface_f0.lammpstrj"
INTERFACE_LOCAL = "input/traj/interface_f0.lammpstrj"
NUCLEATION_REPO = "d-SEAMS/nucleation-frames"
NUCLEATION_REV = ""
NUCLEATION_REMOTE = "data/frames/nucleation_three.lammpstrj"
NUCLEATION_LOCAL = "input/traj/nucleation_three.lammpstrj"


def run_small(root, seams, libdir):
    sys.path.insert(0, str(root / "docs" / "tutorials"))
    import test_brine
    import test_interface
    import test_nucleation

    saved = sys.argv
    try:
        sys.argv = ["test_interface.py", str(root), str(seams), str(libdir)]
        test_interface.main()
        sys.argv = ["test_nucleation.py", str(root), str(libdir)]
        test_nucleation.main()
        sys.argv = ["test_brine.py", str(root), str(libdir)]
        test_brine.main()
    finally:
        sys.argv = saved
    print("small frames ok", PIN_REV)


def fetch_pinned(root, repo, rev, remote, local):
    """Return a local path. A missing revision or a failed fetch uses ``local``."""
    shipped = root / local
    if not shipped.is_file():
        raise SystemExit("missing shipped frame {}".format(shipped))
    if not rev:
        print("{} revision is not pinned; using {}".format(repo, local))
        return shipped
    dest = root / "docs" / ".cache" / repo.split("/", 1)[1] / rev / Path(remote).name
    url = "https://huggingface.co/datasets/{}/resolve/{}/{}".format(repo, rev, remote)
    dest.parent.mkdir(parents=True, exist_ok=True)
    try:
        urllib.request.urlretrieve(url, dest)
    except (urllib.error.URLError, TimeoutError, OSError) as exc:
        print("{} fetch failed ({}); using {}".format(repo, exc, local))
        return shipped
    if dest.stat().st_size < 32:
        print("{} fetch was empty; using {}".format(repo, local))
        return shipped
    print("{} pin {}".format(repo, dest))
    return dest


def fetch_brine(root):
    dest = root / "docs" / ".cache" / "brine-front-230K" / PIN_REV / Path(PIN_FILE).name
    if dest.is_file() and hashlib.sha256(dest.read_bytes()).hexdigest() == PIN_SHA256:
        print("brine pin", dest)
        return dest
    dest.parent.mkdir(parents=True, exist_ok=True)
    try:
        urllib.request.urlretrieve(PIN_URL, dest)
    except (urllib.error.URLError, TimeoutError, OSError) as exc:
        raise SystemExit("could not fetch {}: {}".format(PIN_URL, exc))
    got = hashlib.sha256(dest.read_bytes()).hexdigest()
    if got != PIN_SHA256:
        dest.unlink(missing_ok=True)
        raise SystemExit("brine pin checksum {}".format(got))
    print("brine pin", dest)
    return dest


def main(argv):
    small = "--small" in argv
    args = [arg for arg in argv if arg != "--small"]
    if len(args) != 3:
        raise SystemExit("usage: test_hf_frames.py [--small] ROOT SEAMS LIBDIR")
    root, seams, libdir = Path(args[0]), Path(args[1]), Path(args[2])
    run_small(root, seams, libdir)
    if small:
        print("hf frame tutorials small mode")
        return 0
    for repo, rev, remote, local in (
        (INTERFACE_REPO, INTERFACE_REV, INTERFACE_REMOTE, INTERFACE_LOCAL),
        (NUCLEATION_REPO, NUCLEATION_REV, NUCLEATION_REMOTE, NUCLEATION_LOCAL),
    ):
        fetched = fetch_pinned(root, repo, rev, remote, local)
        shipped = root / local
        if fetched.resolve() != shipped.resolve():
            print("{} fetched {}; the checked frame stays {}".format(repo, fetched, local))
    fetch_brine(root)
    print("hf frame tutorials full mode")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
