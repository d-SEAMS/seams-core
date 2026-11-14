"""Airspeed Velocity arms for PLUMED Q6 against d-SEAMS q6 in PLUMED.

The median over eight samples and the per-atom agreement are the
Snakemake rule ``h2h_plumed_q6``, which runs ``scripts/h2h_plumed_q6.py``.
These arms reuse its frames: set H2H_PLUMED_WORK to the directory given
to that script as --work, and SEAMS_PLUMED_PLUGIN to libSeamsPlumed.so
(default build/engines/libSeamsPlumed.so). Without plumed, the plugin or
the frames, the arms are skipped.

track_s_per_frame is PLUMED's DETAILED_TIMERS for the method's actions,
per frame, so the driver's file reading is not in it.
"""

import importlib.util
import os
import shutil
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
# Loaded by path: this module has the script's name.
_spec = importlib.util.spec_from_file_location("h2h_plumed_q6_script",
                                               ROOT / "scripts" / "h2h_plumed_q6.py")
h2h = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(h2h)

CELLS = {4096: 8, 110592: 24}


class PlumedQ6:
    params = [list(h2h.METHODS), [1, 8], [4096, 110592]]
    param_names = ["method", "threads", "molecules"]
    timeout = 1200

    def setup(self, method, threads, molecules):
        plugin = Path(os.environ.get("SEAMS_PLUMED_PLUGIN",
                                     ROOT / "build" / "engines" / "libSeamsPlumed.so"))
        frames = Path(os.environ.get("H2H_PLUMED_WORK", "")) / f"cells{CELLS[molecules]}" / "frames.xyz"
        if shutil.which("plumed") is None or not plugin.is_file() or not frames.is_file():
            raise NotImplementedError("plumed, libSeamsPlumed.so or the frames are missing")
        self.plugin = plugin
        self.frames = frames
        self.nframes = sum(1 for line in frames.open() if line.startswith("X ")) // molecules
        lines, self.labels, arg = h2h.METHODS[method]
        sw = h2h.SWITCH.format(c=3.5)
        self.deck = [line.format(n=molecules, sw=sw, c=3.5) for line in lines]
        self.deck.append(f"PRINT ARG={arg} FILE=colvar")
        self.env = dict(os.environ, PLUMED_NUM_THREADS=str(threads), OMP_NUM_THREADS=str(threads))

    def track_s_per_frame(self, method, threads, molecules):
        with tempfile.TemporaryDirectory(prefix="asv-plumed-q6-") as tmp:
            log = h2h.driver(self.plugin, self.deck, Path(tmp), self.env, self.frames, "time")
        return h2h.action_seconds(log, self.labels) / self.nframes

    track_s_per_frame.unit = "seconds"
