from pydseamslib._core import PointDouble, PointCloudDouble
import numpy as np
import ase


def map_LAMMPS_IDs_to_atomic_symbols(dict_map: dict(), atms_a: ase.Atoms):
    for pt in atms_a:
        if pt.number in dict_map.keys():
            pt.symbol = dict_map[pt.number]
    return atms_a


def swap_key_value(my_dict_a):
    new_dict = {value: key for key, value in my_dict_a.items()}
    return new_dict


def list_pts(atms: ase.Atoms, lmp_a: dict, splicemask: list[bool], molids: list[int]):
    ptslist = []
    inslice_atms = atms[splicemask]
    map_lmp = swap_key_value(lmp_a)
    all_atm_id = np.asarray([x.index + 1 for x in atms])[splicemask]
    for idx, atm in enumerate(inslice_atms):
        pd = PointDouble()
        pd.x = atm.position[0]
        pd.y = atm.position[1]
        pd.z = atm.position[2]
        pd.c_type = map_lmp[atm.symbol]
        pd.inSlice = True
        pd.atomID = all_atm_id[idx]
        pd.molID = molids[idx]
        ptslist.append(pd)
    return ptslist


def to_pointcloud(
    atms_a: ase.Atoms,
    lmp_a: dict,
    splicemask: list[bool],
    molids: list[int],
    fname,
    currentFrame=[1],
):
    _pcd = PointCloudDouble()
    inslice_atms = atms_a[splicemask]
    _pcd.box = np.diag(np.asarray(inslice_atms.get_cell()))
    _pcd.nop = len(inslice_atms)
    _pcd.pts = list_pts(atms_a, lmp_a, splicemask, molids)
    _pcd.idIndexMap = make_indexmap(atms_a, splicemask)
    _pcd.currentFrame = 1
    _pcd.boxLow = boxlow(fname)
    return _pcd


def make_indexmap(atms: ase.Atoms, splicemask: list[bool]):
    inslice_atms = atms[splicemask]
    atm_id = np.asarray([x.index + 1 for x in atms])[splicemask]
    ret = dict()
    for idx in range(len(inslice_atms)):
        key = int(atm_id[idx])
        ret[key] = idx
    return ret


def boxlow(fname):
    opfile = open(fname)
    data = opfile.read()
    filelist = []
    for i in range(5, 8):
        parfile = data.splitlines()[i]
        needfile = parfile.split()[0]
        floatfile = float(needfile)
        filelist.append(floatfile)
    opfile.close()
    return filelist
