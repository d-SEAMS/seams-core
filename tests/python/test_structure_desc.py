"""Structure descriptor tests on synthetic FCC and BCC lattices.

Mirrors tests/test_structure_desc.cpp: IRA/Horn template assignment,
SOAP spectra, and a linear classifier on Voronoi q4/q6/q8.
"""

from math import isfinite, sqrt

from pydseamslib import _core


def _lattice_cloud(basis, reps, lattice):
    """Build a periodic lattice point cloud of `reps` unit cells per side."""
    cloud = _core.PointCloudDouble()
    pts = []
    id_index_map = {}
    atom_id = 1
    for i in range(reps):
        for j in range(reps):
            for k in range(reps):
                for bx, by, bz in basis:
                    pt = _core.PointDouble()
                    pt.c_type = 1
                    pt.atomID = atom_id
                    pt.molID = atom_id
                    pt.x = (i + bx) * lattice
                    pt.y = (j + by) * lattice
                    pt.z = (k + bz) * lattice
                    pts.append(pt)
                    id_index_map[atom_id] = atom_id - 1
                    atom_id += 1

    cloud.pts = pts
    cloud.nop = len(pts)
    cloud.currentFrame = 1
    box_length = reps * lattice
    cloud.box = [box_length, box_length, box_length]
    cloud.boxLow = [0.0, 0.0, 0.0]
    cloud.idIndexMap = id_index_map
    return cloud


def _fcc_cloud(reps=3, lattice=4.0):
    return _lattice_cloud(
        [
            (0.0, 0.0, 0.0),
            (0.5, 0.5, 0.0),
            (0.5, 0.0, 0.5),
            (0.0, 0.5, 0.5),
        ],
        reps,
        lattice,
    )


def _bcc_cloud(reps=3, lattice=4.0):
    return _lattice_cloud(
        [
            (0.0, 0.0, 0.0),
            (0.5, 0.5, 0.5),
        ],
        reps,
        lattice,
    )


def test_templates_assign_fcc_and_bcc():
    fcc = _fcc_cloud()
    bcc = _bcc_cloud()
    fcc_n = _core.neighListO(rcutoff=3.2, yCloud=fcc, typeI=1)
    bcc_n = _core.neighListO(rcutoff=4.0, yCloud=bcc, typeI=1)

    fcc_hits = _core.classifyTemplates(fcc, fcc_n, 12)
    bcc_hits = _core.classifyTemplates(bcc, bcc_n, 8)

    fcc_ok = sum(1 for h in fcc_hits if h.name == "fcc" and h.rmsd < 0.2)
    bcc_ok = sum(1 for h in bcc_hits if h.name == "bcc" and h.rmsd < 0.2)

    assert fcc_ok > fcc.nop / 2
    assert bcc_ok > bcc.nop / 2


def test_soap_spectrum_size_and_cosine():
    cloud = _fcc_cloud()
    n_list = _core.neighListO(rcutoff=3.2, yCloud=cloud, typeI=1)
    a = list(_core.soapSpectrum(cloud, 0, n_list, 3, 6, 3.2))
    b = list(_core.soapSpectrum(cloud, 1, n_list, 3, 6, 3.2))

    assert len(a) == 3 * 3 * 7
    assert all(isfinite(x) for x in a)
    na = sum(x * x for x in a)
    nb = sum(x * x for x in b)
    dot = sum(x * y for x, y in zip(a, b))
    assert na > 0.0
    assert dot / sqrt(na * nb) > 0.99


def test_linear_classifier_separates_fcc_bcc_voronoi():
    fcc = _fcc_cloud()
    bcc = _bcc_cloud()

    def pack(cloud, cut):
        q4 = _core.steinhardtQlVoronoi(cloud, cut, 4)
        q6 = _core.steinhardtQlVoronoi(cloud, cut, 6)
        q8 = _core.steinhardtQlVoronoi(cloud, cut, 8)
        return [[q4.ql[i], q6.ql[i], q8.ql[i]] for i in range(cloud.nop)]

    fcc_x = pack(fcc, 4.8)
    bcc_x = pack(bcc, 5.6)
    X = fcc_x + bcc_x
    y = [0] * len(fcc_x) + [1] * len(bcc_x)

    clf = _core.LinearClassifier()
    clf.labels = ["fcc", "bcc"]
    clf.fit(X, y)

    fcc_right = sum(1 for row in fcc_x if clf.predict(row) == 0)
    bcc_right = sum(1 for row in bcc_x if clf.predict(row) == 1)
    assert fcc_right == fcc.nop
    assert bcc_right == bcc.nop
