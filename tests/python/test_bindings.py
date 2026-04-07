"""Basic sanity checks for pydseamslib bindings."""

import pydseamslib
from pydseamslib import _core


def test_version():
    assert hasattr(pydseamslib, "__version__")
    assert pydseamslib.__version__ == "2.0.0"


def test_import_core():
    assert hasattr(_core, "PointDouble")
    assert hasattr(_core, "PointCloudDouble")


def test_point_double_construction():
    pt = _core.PointDouble()
    pt.x = 1.0
    pt.y = 2.0
    pt.z = 3.0
    assert pt.x == 1.0
    assert pt.y == 2.0
    assert pt.z == 3.0


def test_pointcloud_double_construction():
    pcd = _core.PointCloudDouble()
    pcd.nop = 0
    assert pcd.nop == 0


def test_readlammps_exists():
    assert hasattr(_core, "readLammpsTrjreduced")


def test_neighlist_exists():
    assert hasattr(_core, "neighListO")


def test_populate_hbonds_exists():
    assert hasattr(_core, "populateHbonds")


def test_ring_network_exists():
    assert hasattr(_core, "ringNetwork")
