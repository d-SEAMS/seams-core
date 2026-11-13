#!/usr/bin/env python3
"""Every lammpstrj reader recovers the genice sH cell edge and tilts.

The dump bound span on x is 31.05 angstrom. The cell edge is
24.8424228 angstrom, with xy = -6.2106057, xz = 0, yz = 0.
"""

import ctypes
import os
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
sys.path.insert(0, str(ROOT / "docs" / "examples"))
sys.path.insert(0, str(ROOT / "docs" / "scripts"))

from lammps_box import read_dump_cell
import render_scenes
import science_render
import core_counts
import render_molecules
import minimal_py

DUMP = ROOT / "input" / "traj" / "genice_sH.lammpstrj"
LX, LY, LZ = 24.8424228, 10.7570846, 20.1312910
XY = -6.2106057


def close(got, want, label):
    if abs(float(got) - want) > 1e-4:
        raise SystemExit(f"{label} is {got}, expected {want}")


def check_box(origin, edges, tilts, label):
    close(origin[0], 0.0, f"{label} ox")
    close(edges[0], LX, f"{label} lx")
    close(edges[1], LY, f"{label} ly")
    close(edges[2], LZ, f"{label} lz")
    close(tilts[0], XY, f"{label} xy")
    close(tilts[1], 0.0, f"{label} xz")
    close(tilts[2], 0.0, f"{label} yz")
    print(f"{label} lx {edges[0]:.7f} xy {tilts[0]:.7f}")


def check_header():
    source = r"""
#include "seams_lammps_box.h"
#include <stdio.h>
int main(void) {
  FILE *in = fopen(DUMP, "r");
  char line[256];
  double lo[3], hi[3], tilt[3], box9[9];
  int axis = 0, reading = 0;
  if (in == NULL) return 1;
  while (fgets(line, sizeof line, in) != NULL) {
    if (!reading) {
      if (strncmp(line, "ITEM: BOX BOUNDS", 16) == 0) reading = 1;
      continue;
    }
    tilt[axis] = 0.0;
    if (sscanf(line, "%lf %lf %lf", &lo[axis], &hi[axis], &tilt[axis]) < 2) return 1;
    if (++axis == 3) break;
  }
  fclose(in);
  seams_lammps_bounds(lo, hi, tilt, box9);
  printf("%.10f %.10f %.10f %.10f %.10f %.10f %.10f\n",
         box9[0], box9[3], box9[4], box9[5], box9[6], box9[7], box9[8]);
  return 0;
}
""".replace("DUMP", f'"{DUMP}"')
    # strncmp needs string.h
    source = "#include <string.h>\n" + source
    with tempfile.TemporaryDirectory() as tmp:
        path = Path(tmp) / "shcell.c"
        path.write_text(source)
        binary = Path(tmp) / "shcell"
        subprocess.check_call(
            ["/usr/bin/cc", "-I", str(ROOT / "include"), str(path), "-o", str(binary)],
        )
        text = subprocess.check_output([str(binary)], text=True).split()
    origin = (float(text[0]), 0.0, 0.0)
    edges = tuple(float(item) for item in text[1:4])
    tilts = tuple(float(item) for item in text[4:7])
    check_box(origin, edges, tilts, "C header")


def check_lang():
    libdir = os.environ.get("SEAMS_LIBDIR")
    if not libdir:
        candidates = [ROOT / "bbdir" / "bindings"]
        candidates.extend(Path(item) for item in os.environ.get("LD_LIBRARY_PATH", "").split(":") if item)
        for candidate in candidates:
            if (candidate / "libseams_lang.so").is_file():
                libdir = str(candidate)
                break
    if not libdir:
        raise SystemExit("libseams_lang.so is required for the language reader")
    lib = ctypes.CDLL(str(Path(libdir) / "libseams_lang.so"))
    lib.seams_frame_read.argtypes = [ctypes.c_char_p, ctypes.c_int, ctypes.c_int]
    lib.seams_frame_read.restype = ctypes.c_void_p
    lib.seams_frame_nop.argtypes = [ctypes.c_void_p]
    lib.seams_frame_nop.restype = ctypes.c_int
    lib.seams_frame_box.argtypes = [ctypes.c_void_p, ctypes.POINTER(ctypes.c_double)]
    lib.seams_frame_free.argtypes = [ctypes.c_void_p]
    frame = lib.seams_frame_read(str(DUMP).encode(), 1, 1)
    if not frame or lib.seams_frame_nop(frame) != 136:
        raise SystemExit("language reader did not return the sH frame")
    box = (ctypes.c_double * 9)()
    lib.seams_frame_box(frame, box)
    lib.seams_frame_free(frame)
    check_box(box[0:3], box[3:6], box[6:9], "seams_frame_box")


def main():
    origin, edges, tilts = read_dump_cell(DUMP)
    check_box(origin, edges, tilts, "lammps_box")
    vectors = render_scenes._box_vectors(DUMP)
    check_box((0.0, 0.0, 0.0), (vectors[0][0], vectors[1][1], vectors[2][2]),
              (vectors[1][0], vectors[2][0], vectors[2][1]), "render_scenes")
    lengths, _coords = science_render.read_dump_frame(DUMP, 1)
    close(lengths[0], LX, "science_render lx")
    close(lengths[1], LY, "science_render ly")
    close(lengths[2], LZ, "science_render lz")
    print(f"science_render lx {lengths[0]:.7f}")
    _xyz, _types, box = core_counts._load_dump(DUMP)
    check_box(box[0:3], box[3:6], box[6:9], "core_counts")
    frame = render_molecules.read_dump(DUMP)
    check_box(frame.box9[0:3], frame.box9[3:6], frame.box9[6:9], "render_molecules")
    _points, edges = minimal_py.read_dump(DUMP)
    close(edges[0], LX, "minimal_py lx")
    close(edges[1], LY, "minimal_py ly")
    close(edges[2], LZ, "minimal_py lz")
    print(f"minimal_py lx {edges[0]:.7f}")
    check_header()
    check_lang()
    print("sH cell ok")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
