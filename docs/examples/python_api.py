#!/usr/bin/env python3
"""Eight-atom CHILL+ and q4/q6 through dseams, or through libseams."""

import ctypes
import math
import os
import sys
from pathlib import Path


def diamond():
    a = 2.75 * 4.0 / math.sqrt(3.0)
    fcc = ((0, 0, 0), (0.5, 0.5, 0), (0.5, 0, 0.5), (0, 0.5, 0.5))
    basis = ((0, 0, 0), (0.25, 0.25, 0.25))
    xyz = []
    for shift in fcc:
        for site in basis:
            xyz.extend((shift[i] + site[i]) * a for i in range(3))
    return xyz, a


def via_dseams(xyz, edge):
    import numpy as np

    import dseams

    pos = np.asarray(xyz, dtype=np.float64).reshape(8, 3)
    box = np.asarray([edge, edge, edge], dtype=np.float64)
    dseams.set_num_threads(1)
    if dseams.num_threads() != 1:
        raise SystemExit("dseams thread width")
    labels, counts = dseams.chill_plus(pos, box)
    q4, q6 = dseams.q4q6(pos, box)
    if int(counts[0]) != 8 or len(labels) != 8 or len(q6) != 8:
        raise SystemExit("dseams counts")
    if float(np.mean(q6)) <= 0.3:
        raise SystemExit("dseams q6")
    print("python api ok dseams", dseams.version())


def via_ctypes(libdir, xyz, edge):
    lib = ctypes.CDLL(os.path.join(libdir, "libseams.so"))
    lib.seams_version_string.restype = ctypes.c_char_p
    lib.seams_set_num_threads.argtypes = [ctypes.c_int]
    lib.seams_set_num_threads.restype = ctypes.c_int
    lib.seams_get_num_threads.restype = ctypes.c_int
    lib.seams_system_create.restype = ctypes.c_void_p
    lib.seams_system_set_ortho.argtypes = [
        ctypes.c_void_p,
        ctypes.POINTER(ctypes.c_double),
        ctypes.c_int,
        ctypes.c_double,
        ctypes.c_double,
        ctypes.c_double,
        ctypes.c_double,
        ctypes.c_double,
        ctypes.c_double,
        ctypes.c_void_p,
        ctypes.c_void_p,
    ]
    lib.seams_system_set_ortho.restype = ctypes.c_int
    lib.seams_chill_plus.argtypes = [
        ctypes.c_void_p,
        ctypes.c_double,
        ctypes.c_int,
        ctypes.c_int,
        ctypes.c_int,
        ctypes.POINTER(ctypes.c_int),
        ctypes.POINTER(ctypes.c_int),
    ]
    lib.seams_chill_plus.restype = ctypes.c_int
    lib.seams_q4q6.argtypes = [
        ctypes.c_void_p,
        ctypes.c_double,
        ctypes.c_int,
        ctypes.POINTER(ctypes.c_double),
        ctypes.POINTER(ctypes.c_double),
    ]
    lib.seams_q4q6.restype = ctypes.c_int
    lib.seams_system_destroy.argtypes = [ctypes.c_void_p]
    lib.seams_system_count.argtypes = [ctypes.c_void_p]
    lib.seams_system_count.restype = ctypes.c_int
    lib.seams_handle_error.argtypes = [ctypes.c_void_p]
    lib.seams_handle_error.restype = ctypes.c_char_p
    lib.seams_system_read.argtypes = [ctypes.c_void_p, ctypes.c_char_p, ctypes.c_int, ctypes.c_int]
    lib.seams_system_read.restype = ctypes.c_int
    lib.seams_chill.argtypes = [
        ctypes.c_void_p, ctypes.c_double, ctypes.c_int, ctypes.c_int, ctypes.c_int,
        ctypes.POINTER(ctypes.c_int), ctypes.c_size_t, ctypes.POINTER(ctypes.c_int),
    ]
    lib.seams_chill.restype = ctypes.c_int
    lib.seams_steinhardt.argtypes = [
        ctypes.c_void_p, ctypes.c_double, ctypes.c_int, ctypes.c_int,
        ctypes.POINTER(ctypes.c_double), ctypes.c_size_t,
        ctypes.POINTER(ctypes.c_double), ctypes.c_size_t,
    ]
    lib.seams_steinhardt.restype = ctypes.c_int
    lib.seams_cages.argtypes = [
        ctypes.c_void_p, ctypes.c_double, ctypes.c_int, ctypes.c_int, ctypes.c_int,
        ctypes.c_int, ctypes.POINTER(ctypes.c_void_p),
    ]
    lib.seams_cages.restype = ctypes.c_int
    lib.seams_cages_counts.argtypes = [ctypes.c_void_p] + [ctypes.POINTER(ctypes.c_int)] * 4
    lib.seams_cages_counts.restype = ctypes.c_int
    lib.seams_cages_destroy.argtypes = [ctypes.c_void_p]
    lib.seams_signature_cages.argtypes = [
        ctypes.c_void_p, ctypes.c_double, ctypes.c_int, ctypes.c_int, ctypes.c_char_p,
        ctypes.POINTER(ctypes.c_size_t), ctypes.POINTER(ctypes.c_size_t),
    ]
    lib.seams_signature_cages.restype = ctypes.c_int
    lib.seams_neighbours_knn.argtypes = [
        ctypes.c_void_p, ctypes.c_int, ctypes.c_double, ctypes.c_int, ctypes.c_int,
        ctypes.POINTER(ctypes.c_void_p),
    ]
    lib.seams_neighbours_knn.restype = ctypes.c_int
    lib.seams_neighbours_destroy.argtypes = [ctypes.c_void_p]
    lib.seams_rdf.restype = ctypes.c_int
    lib.seams_running_cn.restype = ctypes.c_int
    lib.seams_coordination.argtypes = [
        ctypes.c_void_p, ctypes.c_double, ctypes.c_int, ctypes.c_int, ctypes.c_int,
        ctypes.POINTER(ctypes.c_double),
    ]
    lib.seams_coordination.restype = ctypes.c_int
    lib.seams_hbonds.argtypes = [
        ctypes.c_void_p, ctypes.c_double, ctypes.c_int, ctypes.c_int, ctypes.c_double,
        ctypes.c_double, ctypes.POINTER(ctypes.c_int),
    ]
    lib.seams_hbonds.restype = ctypes.c_int
    lib.seams_density_axis.restype = ctypes.c_int
    lib.seams_fingerprint.argtypes = [
        ctypes.c_void_p, ctypes.c_double, ctypes.c_int, ctypes.c_int, ctypes.c_int,
        ctypes.c_int, ctypes.c_int, ctypes.POINTER(ctypes.c_void_p),
    ]
    lib.seams_fingerprint.restype = ctypes.c_int
    lib.seams_fingerprint_destroy.argtypes = [ctypes.c_void_p]
    lib.seams_fingerprint_key.argtypes = [
        ctypes.c_void_p, ctypes.c_void_p, ctypes.c_size_t, ctypes.POINTER(ctypes.c_size_t),
    ]
    lib.seams_fingerprint_key.restype = ctypes.c_int
    lib.seams_contact_pairs.argtypes = [
        ctypes.c_void_p, ctypes.c_char_p, ctypes.POINTER(ctypes.c_size_t),
        ctypes.POINTER(ctypes.c_int), ctypes.POINTER(ctypes.c_int),
    ]
    lib.seams_contact_pairs.restype = ctypes.c_int
    lib.seams_domains.argtypes = [
        ctypes.c_void_p, ctypes.c_double, ctypes.c_char_p, ctypes.c_int,
        ctypes.POINTER(ctypes.c_int), ctypes.POINTER(ctypes.c_int), ctypes.POINTER(ctypes.c_double),
    ]
    lib.seams_domains.restype = ctypes.c_int
    lib.seams_prisms.argtypes = [
        ctypes.c_void_p, ctypes.c_double, ctypes.c_int, ctypes.POINTER(ctypes.c_int),
        ctypes.POINTER(ctypes.c_int), ctypes.c_size_t,
    ]
    lib.seams_prisms.restype = ctypes.c_int
    lib.seams_monolayer.argtypes = [
        ctypes.c_void_p, ctypes.c_double, ctypes.c_int, ctypes.c_int, ctypes.c_double,
        ctypes.POINTER(ctypes.c_int), ctypes.c_size_t,
    ]
    lib.seams_monolayer.restype = ctypes.c_int
    lib.seams_rings_find.argtypes = [
        ctypes.c_void_p, ctypes.c_double, ctypes.c_int, ctypes.c_int, ctypes.c_int,
        ctypes.POINTER(ctypes.c_void_p),
    ]
    lib.seams_rings_find.restype = ctypes.c_int
    lib.seams_rings_destroy.argtypes = [ctypes.c_void_p]
    lib.seams_rings_vertices.argtypes = [
        ctypes.c_void_p, ctypes.c_int, ctypes.POINTER(ctypes.c_int), ctypes.c_size_t,
        ctypes.POINTER(ctypes.c_size_t),
    ]
    lib.seams_rings_vertices.restype = ctypes.c_int
    lib.seams_site_table.argtypes = [
        ctypes.c_void_p, ctypes.c_char_p, ctypes.POINTER(ctypes.c_int), ctypes.c_size_t,
    ]
    lib.seams_site_table.restype = ctypes.c_int
    if lib.seams_set_num_threads(1) != 0 or lib.seams_get_num_threads() != 1:
        raise SystemExit("ctypes thread width")
    points = (ctypes.c_double * len(xyz))(*xyz)
    labels = (ctypes.c_int * 8)()
    counts = (ctypes.c_int * 9)()
    q4 = (ctypes.c_double * 8)()
    q6 = (ctypes.c_double * 8)()
    sys_handle = lib.seams_system_create()
    if lib.seams_system_set_ortho(sys_handle, points, 8, 0, 0, 0, edge, edge, edge, None, None) != 0:
        raise SystemExit("ctypes frame")
    if lib.seams_chill_plus(sys_handle, 3.5, 0, 0, 0, labels, counts) != 0:
        raise SystemExit("ctypes chill+")
    if lib.seams_q4q6(sys_handle, 3.5, 0, q4, q6) != 0:
        raise SystemExit("ctypes q4q6")
    mean = sum(q6) / 8.0
    if counts[0] != 8 or mean <= 0.3:
        raise SystemExit("ctypes result")
    exercise(lib, sys_handle)
    lib.seams_system_destroy(sys_handle)
    root = Path(__file__).resolve().parents[2]
    expect = {
        "input/traj/ice1h.lammpstrj": (128, -1),
        "input/traj/genice_sI.lammpstrj": (46, 1),
        "input/xyz/configuration.xyz": (256, -1),
        "input/con/tiny_multi_cuh2.con": (None, -1),
        "input/traj/genice_sH.lammpstrj": (136, 1),
    }
    for rel, (count, filt) in expect.items():
        got = read_count(lib, str(root / rel), filt)
        if count is not None and got != count:
            raise SystemExit("ctypes read {} {}".format(rel, got))
        if got <= 0:
            raise SystemExit("ctypes read empty {}".format(rel))
    guest = lib.seams_system_create()
    guest_path = str(root / "input/traj/genice_sI.lammpstrj").encode()
    if lib.seams_system_read(guest, guest_path, 1, 1) != 0:
        fail(lib, guest, "sI")
    n_cages = ctypes.c_size_t()
    n_atoms = ctypes.c_size_t()
    if lib.seams_signature_cages(guest, 3.5, 0, 0, b"512", n_cages, n_atoms) != 0:
        fail(lib, guest, "signature")
    if n_cages.value != 2 or n_atoms.value != 40:
        raise SystemExit("ctypes signature {} {}".format(n_cages.value, n_atoms.value))
    lib.seams_system_destroy(guest)
    print("python api ok ctypes", lib.seams_version_string().decode())


def fail(lib, handle, label):
    text = lib.seams_handle_error(handle)
    detail = text.decode() if text else ""
    raise SystemExit("ctypes {} {}".format(label, detail))


def exercise(lib, sys_handle):
    labels = (ctypes.c_int * 8)()
    counts = (ctypes.c_int * 9)()
    if lib.seams_chill(sys_handle, 3.5, 0, 0, 0, labels, 8, counts) != 0:
        fail(lib, sys_handle, "chill")
    ql = (ctypes.c_double * 8)()
    bar = (ctypes.c_double * 8)()
    if lib.seams_steinhardt(sys_handle, 3.5, 0, 6, ql, 8, bar, 8) != 0:
        fail(lib, sys_handle, "steinhardt")
    if sum(bar) / 8.0 <= 0.0:
        raise SystemExit("ctypes lechner-dellago")
    cages = ctypes.c_void_p()
    if lib.seams_cages(sys_handle, 3.5, 0, 2, 4, 0, ctypes.byref(cages)) != 0:
        fail(lib, sys_handle, "cages")
    n_hex = ctypes.c_int()
    n_cubic = ctypes.c_int()
    n_water = ctypes.c_int()
    n_both = ctypes.c_int()
    if lib.seams_cages_counts(cages, n_hex, n_cubic, n_water, n_both) != 0:
        fail(lib, cages, "cage counts")
    lib.seams_cages_destroy(cages)
    knn = ctypes.c_void_p()
    if lib.seams_neighbours_knn(sys_handle, 4, 3.5, 0, 1, ctypes.byref(knn)) != 0:
        fail(lib, sys_handle, "knn")
    lib.seams_neighbours_destroy(knn)
    cn = ctypes.c_double()
    if lib.seams_coordination(sys_handle, 3.5, 8, 0, 0, cn) != 0:
        fail(lib, sys_handle, "coordination")
    n_bonds = ctypes.c_int()
    if lib.seams_hbonds(sys_handle, 3.5, 1, 2, 2.42, 30.0, n_bonds) != 0:
        fail(lib, sys_handle, "hbonds")
    n_out = ctypes.c_size_t()
    r = (ctypes.c_double * 16)()
    g = (ctypes.c_double * 16)()
    lib.seams_rdf.argtypes = [
        ctypes.c_void_p, ctypes.c_double, ctypes.c_int, ctypes.c_int, ctypes.c_int,
        ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double), ctypes.c_void_p,
        ctypes.c_size_t, ctypes.POINTER(ctypes.c_size_t), ctypes.POINTER(ctypes.c_double),
    ]
    volume = ctypes.c_double()
    if lib.seams_rdf(sys_handle, 6.0, 8, 0, 0, r, g, None, 16, n_out, volume) != 0:
        fail(lib, sys_handle, "rdf")
    lib.seams_running_cn.argtypes = [
        ctypes.c_void_p, ctypes.c_double, ctypes.c_int, ctypes.c_int, ctypes.c_int,
        ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double), ctypes.c_size_t,
        ctypes.POINTER(ctypes.c_size_t),
    ]
    if lib.seams_running_cn(sys_handle, 6.0, 8, 0, 0, r, g, 16, n_out) != 0:
        fail(lib, sys_handle, "running_cn")
    z = (ctypes.c_double * 8)()
    rho = (ctypes.c_double * 8)()
    lib.seams_density_axis.argtypes = [
        ctypes.c_void_p, ctypes.c_int, ctypes.c_int, ctypes.c_int,
        ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double), ctypes.c_size_t,
        ctypes.POINTER(ctypes.c_size_t),
    ]
    if lib.seams_density_axis(sys_handle, 0, 4, 2, z, rho, 8, n_out) != 0:
        fail(lib, sys_handle, "density")
    fp = ctypes.c_void_p()
    if lib.seams_fingerprint(sys_handle, 3.5, 0, 4, 0, 2, 0, ctypes.byref(fp)) != 0:
        fail(lib, sys_handle, "fingerprint")
    if lib.seams_fingerprint_key(fp, None, 0, n_out) != 0 or n_out.value <= 0:
        fail(lib, fp, "fingerprint key")
    lib.seams_fingerprint_destroy(fp)
    n_pairs = ctypes.c_size_t()
    n_cation = ctypes.c_int()
    n_anion = ctypes.c_int()
    if lib.seams_contact_pairs(sys_handle, b"1=cationHead", n_pairs, n_cation, n_anion) != 0:
        fail(lib, sys_handle, "pairs")
    n_subset = ctypes.c_int()
    largest = ctypes.c_int()
    percolation = ctypes.c_double()
    if lib.seams_domains(sys_handle, 3.5, b"1=waterO", 0, n_subset, largest, percolation) != 0:
        fail(lib, sys_handle, "domains")
    perfect = (ctypes.c_int * 7)()
    imperfect = (ctypes.c_int * 7)()
    if lib.seams_prisms(sys_handle, 3.5, 0, perfect, imperfect, 7) != 0:
        fail(lib, sys_handle, "prisms")
    sheet = (ctypes.c_int * 7)()
    if lib.seams_monolayer(sys_handle, 3.5, 0, 6, 100.0, sheet, 7) != 0:
        fail(lib, sys_handle, "monolayer")
    rings = ctypes.c_void_p()
    if lib.seams_rings_find(sys_handle, 3.5, 0, 0, 6, ctypes.byref(rings)) != 0:
        fail(lib, sys_handle, "rings")
    atoms = (ctypes.c_int * 12)()
    if lib.seams_rings_vertices(rings, 0, atoms, 12, n_out) != 0 or n_out.value <= 0:
        fail(lib, rings, "vertices")
    lib.seams_rings_destroy(rings)
    kinds = (ctypes.c_int * 8)()
    if lib.seams_site_table(sys_handle, b"1=waterO", kinds, 8) != 0:
        fail(lib, sys_handle, "site")


def read_count(lib, path, type_filter):
    handle = lib.seams_system_create()
    code = lib.seams_system_read(handle, path.encode(), 1, type_filter)
    if code != 0:
        fail(lib, handle, path)
    count = lib.seams_system_count(handle)
    lib.seams_system_destroy(handle)
    return count


def main():
    xyz, edge = diamond()
    loaded = False
    try:
        via_dseams(xyz, edge)
        loaded = True
    except ImportError:
        loaded = False
    if len(sys.argv) >= 2:
        via_ctypes(sys.argv[1], xyz, edge)
    elif not loaded:
        raise SystemExit("pass the directory that holds libseams.so")


if __name__ == "__main__":
    main()
