"""User-facing one-frame handle. Load a file, ASE Atoms, or arrays; ask for ice."""

from __future__ import annotations

from collections import Counter
from pathlib import Path

from . import _core


class IceCounts(dict):
    """CHILL / CHILL+ histogram. Use counts.cubic or counts['cubic']."""

    def __getattr__(self, name):
        if name.startswith("_"):
            raise AttributeError(name)
        return self.get(name, 0)

    def __repr__(self):
        parts = [f"{k}={v}" for k, v in sorted(self.items()) if v]
        return "IceCounts(" + ", ".join(parts) + ")"


class CageScore:
    """Per-atom HC/DDC membership. A molecule in an HC is ice Ih, DDC is Ic."""

    def __init__(self, hc, ddc):
        self.hc = list(hc)
        self.ddc = list(ddc)

    @property
    def n_ih(self):
        return sum(1 for flag in self.hc if flag)

    @property
    def n_ic(self):
        return sum(1 for flag in self.ddc if flag)

    @property
    def n_water(self):
        return sum(
            1 for h, d in zip(self.hc, self.ddc) if not h and not d
        )

    def __repr__(self):
        return (
            f"CageScore(n_ih={self.n_ih}, n_ic={self.n_ic}, "
            f"n_water={self.n_water})"
        )


def _cloud_from_positions(positions, cell, numbers, box_low=None):
    n = len(positions)
    if n == 0:
        raise ValueError("no atoms to load")
    if len(cell) != 3:
        raise ValueError("cell must be three box lengths [lx, ly, lz]")
    cloud = _core.PointCloudDouble()
    cloud.nop = n
    cloud.currentFrame = 1
    cloud.box = [float(cell[0]), float(cell[1]), float(cell[2])]
    if box_low is None:
        cloud.boxLow = [0.0, 0.0, 0.0]
    else:
        cloud.boxLow = [float(box_low[0]), float(box_low[1]), float(box_low[2])]
    pts = []
    id_map = {}
    for i, xyz in enumerate(positions):
        pt = _core.PointDouble()
        pt.x = float(xyz[0])
        pt.y = float(xyz[1])
        pt.z = float(xyz[2])
        pt.c_type = int(numbers[i]) if numbers is not None else 1
        pt.atomID = i + 1
        pt.molID = i + 1
        pt.inSlice = True
        pts.append(pt)
        id_map[i + 1] = i
    cloud.pts = pts
    cloud.idIndexMap = id_map
    return cloud


def _guess_lammps_type(filename, frame, region):
    low, high = region if region is not None else ([0, 0, 0], [0, 0, 0])
    sliced = region is not None
    last = None
    for type_i in (2, 1):
        cloud = _core.readLammpsTrjreduced(
            filename=filename,
            targetFrame=frame,
            typeI=type_i,
            isSlice=sliced,
            coordLow=list(low),
            coordHigh=list(high),
        )
        if cloud.nop > 0:
            return cloud, type_i
    raise ValueError(f"{filename} frame {frame} has no atoms of type 1 or 2")


class Frame:
    """One configuration: neighbours, rings, CHILL(+), and cage membership.

    Load a LAMMPS dump, an ASE Atoms, or raw arrays. Then call chill_plus()
    or cages(). Classification does not write files.

    Parameters
    ----------
    filename : path
        LAMMPS dump. Types 1 and 2 are treated as hydrogen and oxygen
        unless atom_type is set.
    frame : int
        1-indexed frame.
    atom_type : int or None
        Species to analyse. None picks oxygen (type 2) if that type is
        present, otherwise type 1 (mW-style single-site dumps).
    cutoff : float
        Neighbour cutoff in Angstroms.
    bonded : {"auto", "hbond", "cutoff"}
        Graph for rings. "auto" uses hydrogen bonds when hydrogens are
        available, otherwise the cutoff neighbour list.
    region : ((xlo,ylo,zlo), (xhi,yhi,zhi)) or None
    """

    def __init__(
        self,
        filename=None,
        frame=1,
        atom_type=2,
        cutoff=3.5,
        bonded="hbond",
        region=None,
        *,
        cloud=None,
        h_cloud=None,
        symbols=None,
    ):
        if bonded not in ("hbond", "cutoff", "auto"):
            raise ValueError('bonded must be "hbond", "cutoff", or "auto"')
        self.region = region
        self.filename = (
            str(Path(filename).resolve()) if filename is not None else None
        )
        self.frame = frame
        self.cutoff = cutoff
        self._h_cloud = h_cloud
        self._symbols = symbols
        self._nlist = None
        self._hbonds = None
        self._rings = None
        self._classifier = None
        self._ring_updater = _core.RingUpdater(6)
        self._affiliation_updater = None

        if cloud is not None:
            self.cloud = cloud
            self.atom_type = atom_type
        elif self.filename is not None:
            if atom_type is None:
                self.cloud, self.atom_type = _guess_lammps_type(
                    self.filename, frame, region
                )
            else:
                self.atom_type = atom_type
                self.cloud = self._read(frame)
        else:
            raise TypeError("Frame needs a filename, cloud=, from_ase, or from_arrays")

        if bonded == "auto":
            self.bonded = "hbond" if self._can_hbond() else "cutoff"
        else:
            self.bonded = bonded

    def _can_hbond(self):
        if self._h_cloud is not None:
            return self._h_cloud.nop > 0
        return self.filename is not None and self.atom_type == 2

    def _read(self, frame):
        low, high = self.region if self.region is not None else (
            [0, 0, 0], [0, 0, 0]
        )
        return _core.readLammpsTrjreduced(
            filename=self.filename,
            targetFrame=frame,
            typeI=self.atom_type,
            isSlice=self.region is not None,
            coordLow=list(low),
            coordHigh=list(high),
        )

    @classmethod
    def from_file(cls, filename, frame=1, atom_type=None, cutoff=3.5,
                  bonded="auto", region=None):
        """Load a LAMMPS dump. atom_type None means oxygen if present."""
        return cls(
            filename,
            frame=frame,
            atom_type=atom_type if atom_type is not None else 2,
            cutoff=cutoff,
            bonded=bonded,
            region=region,
        ) if atom_type is not None else cls._from_file_guess(
            filename, frame, cutoff, bonded, region
        )

    @classmethod
    def _from_file_guess(cls, filename, frame, cutoff, bonded, region):
        path = str(Path(filename).resolve())
        cloud, type_i = _guess_lammps_type(path, frame, region)
        return cls(
            filename=path,
            frame=frame,
            atom_type=type_i,
            cutoff=cutoff,
            bonded=bonded,
            region=region,
            cloud=cloud,
        )

    @classmethod
    def from_arrays(cls, positions, cell, numbers=None, cutoff=3.5,
                    bonded="cutoff", box_low=None):
        """Build a frame from (N, 3) positions and three box lengths."""
        n = len(positions)
        nums = list(numbers) if numbers is not None else [1] * n
        types = sorted(set(int(x) for x in nums))
        atom_type = types[0] if len(types) == 1 else 1
        cloud = _cloud_from_positions(positions, cell, nums, box_low)
        return cls(
            atom_type=atom_type,
            cutoff=cutoff,
            bonded=bonded,
            cloud=cloud,
        )

    @classmethod
    def from_ase(cls, atoms, select="O", cutoff=3.5, bonded="auto"):
        """Load an ASE Atoms. select is a symbol, atomic number, or None."""
        from .aseio import frame_from_ase

        return frame_from_ase(cls, atoms, select=select, cutoff=cutoff,
                              bonded=bonded)

    def to_ase(self):
        """ASE Atoms for this frame. Adds arrays['ice_type'] after CHILL."""
        from .aseio import frame_to_ase

        return frame_to_ase(self)

    @property
    def n_atoms(self):
        return self.cloud.nop

    @property
    def box(self):
        return list(self.cloud.box)

    @property
    def positions(self):
        return [(pt.x, pt.y, pt.z) for pt in self.cloud.pts]

    @property
    def neighbor_list(self):
        if self._nlist is None:
            self._nlist = _core.neighListO(
                rcutoff=self.cutoff, yCloud=self.cloud, typeI=self.atom_type
            )
        return self._nlist

    @property
    def hbonds(self):
        if self._hbonds is None:
            if self._h_cloud is not None:
                self._hbonds = _core.populateHbondsWithInputClouds(
                    yCloud=self.cloud,
                    hCloud=self._h_cloud,
                    nList=self.neighbor_list,
                )
            elif self.filename is not None:
                self._hbonds = _core.populateHbonds(
                    filename=self.filename,
                    yCloud=self.cloud,
                    nList=self.neighbor_list,
                    targetFrame=self.frame,
                    Htype=1,
                )
            else:
                raise ValueError(
                    "Hydrogen-bond analysis needs hydrogens. Pass ASE Atoms "
                    "that include H, a LAMMPS dump with H, or set bonded='cutoff'."
                )
        return self._hbonds

    def load_frame(self, frame):
        if self.filename is None:
            raise ValueError("load_frame needs a trajectory file")
        self.frame = frame
        self.cloud = self._read(frame)
        self._nlist = None
        self._hbonds = None
        self._rings = None

    @property
    def bonds_by_index(self):
        source = (
            self.neighbor_list if self.bonded == "cutoff" else self.hbonds
        )
        return _core.neighbourListByIndex(yCloud=self.cloud, nList=source)

    @property
    def rings(self):
        if self._rings is None:
            self._rings = self._ring_updater.update(self.bonds_by_index)
        return self._rings

    @property
    def rings_recomputed_sources(self):
        return self._ring_updater.lastRecomputedSources()

    def _count_ice(self):
        names = [pt.iceType.name for pt in self.cloud.pts]
        return IceCounts(Counter(names))

    def chill_plus(self):
        """CHILL+ labels. Does not write a file."""
        _core.getCorrelPlus(
            yCloud=self.cloud, nList=self.neighbor_list, isSlice=False
        )
        _core.getIceTypePlusNoPrint(
            yCloud=self.cloud, nList=self.neighbor_list, isSlice=False
        )
        return self._count_ice()

    def chill(self):
        """CHILL labels. Does not write a file."""
        _core.getCorrel(
            yCloud=self.cloud, nList=self.neighbor_list, isSlice=False
        )
        _core.getIceTypeNoPrint(
            yCloud=self.cloud, nList=self.neighbor_list, isSlice=False
        )
        return self._count_ice()

    def classify_chill_plus(self):
        return self.chill_plus()

    def classify_chill(self):
        return self.chill()

    def cages(self, seeded=True, k=4, candidate_cutoff=None):
        """Ice score: HC = Ih, DDC = Ic, neither = water.

        seeded=True is the hysteresis construction (mutual four-nearest
        seeds, union-graph completion). seeded=False is cutoff-graph
        affiliation on this frame's six-rings.
        """
        if seeded:
            return self.seeded_affiliation(k=k, candidate_cutoff=candidate_cutoff)
        aff = self.cage_affiliation()
        n = self.n_atoms
        hc = [False] * n
        ddc = [False] * n
        for ring, is_hc, is_ddc in zip(aff["six_rings"], aff["hc"], aff["ddc"]):
            for atom in ring:
                if 0 <= atom < n:
                    hc[atom] = hc[atom] or bool(is_hc)
                    ddc[atom] = ddc[atom] or bool(is_ddc)
        return CageScore(hc, ddc)

    def cage_affiliation(self):
        six = [r for r in self.rings if len(r) == 6]
        if self._affiliation_updater is None:
            self._affiliation_updater = _core.AffiliationUpdater()
        hc, ddc = self._affiliation_updater.update(six, self.bonds_by_index)
        return {
            "six_rings": six,
            "hc": list(hc),
            "ddc": list(ddc),
            "reclassified": self._affiliation_updater.lastReclassified(),
        }

    def seeded_affiliation(self, k=4, candidate_cutoff=None):
        cut = (
            self.cutoff + 1.5 if candidate_cutoff is None else candidate_cutoff
        )
        strict = _core.neighbourListByIndex(
            self.cloud,
            _core.kNearestNeighbourList(
                self.cloud, k, cut, self.atom_type, True
            ),
        )
        union = _core.neighbourListByIndex(
            self.cloud,
            _core.kNearestNeighbourList(
                self.cloud, k, cut, self.atom_type, False
            ),
        )
        six_s = [r for r in _core.ringNetwork(strict, 6) if len(r) == 6]
        six_u = [r for r in _core.ringNetwork(union, 6) if len(r) == 6]
        hc, ddc = _core.seededCageAffiliation(six_s, strict, six_u, union)
        return CageScore(hc, ddc)

    def find_prisms(self, output_dir="output/", max_depth=6, shape_matching=False):
        hbonds_idx = _core.neighbourListByIndex(
            yCloud=self.cloud, nList=self.hbonds
        )
        _core.prismAnalysis(
            path=output_dir,
            rings=self.rings,
            nList=hbonds_idx,
            yCloud=self.cloud,
            maxDepth=max_depth,
            atomID=0,
            firstFrame=self.frame,
            currentFrame=self.frame,
            doShapeMatching=shape_matching,
        )

    def monolayer_rings(self, output_dir, sheet_area, max_depth=4):
        rings = _core.ringNetwork(self.bonds_by_index, max_depth)
        _core.polygonRingAnalysis(
            path=str(output_dir) + "/",
            rings=rings,
            nList=self.bonds_by_index,
            yCloud=self.cloud,
            maxDepth=max_depth,
            sheetArea=sheet_area,
            firstFrame=self.frame,
        )
        counts = {}
        cov = (
            Path(output_dir) / "topoMonolayer" / "coverageAreaXY.dat"
        ).read_text().splitlines()
        fields = cov[-1].split()[1:]
        for size, n, area in zip(fields[::3], fields[1::3], fields[2::3]):
            counts[int(size)] = {"count": int(n), "coverage_xy": float(area)}
        return counts

    def rdf_2d(self, output_dir, cutoff=12.0, binwidth=0.05):
        _core.rdf2Danalysis_AA(
            path=str(output_dir) + "/",
            rdfValues=[],
            yCloud=self.cloud,
            cutoff=cutoff,
            binwidth=binwidth,
            firstFrame=self.frame,
            finalFrame=self.frame,
        )
        r, g = [], []
        rdf = (
            Path(output_dir) / "topoMonolayer" / "rdf.dat"
        ).read_text().splitlines()
        for line in rdf:
            parts = line.split()
            if len(parts) == 2:
                r.append(float(parts[0]))
                g.append(float(parts[1]))
        return r, g

    def steinhardt(self, order_l=6):
        result = _core.steinhardtQl(
            yCloud=self.cloud, nList=self.neighbor_list, orderL=order_l
        )
        return {"ql": list(result.ql), "ql_bar": list(result.qlBar)}

    def steinhardt_voronoi(self, order_l=6, cutoff=None):
        result = _core.steinhardtQlVoronoi(
            yCloud=self.cloud,
            candidateCutoff=cutoff if cutoff is not None else self.cutoff,
            orderL=order_l,
        )
        return {"ql": list(result.ql), "ql_bar": list(result.qlBar)}

    def classify_templates(self, k_neigh=12):
        hits = _core.classifyTemplates(self.cloud, self.neighbor_list, k_neigh)
        rows = []
        for h in hits:
            name = h.name
            kind = name or getattr(h.kind, "name", str(h.kind))
            rows.append({"name": name, "rmsd": h.rmsd, "kind": kind})
        return rows

    def soap(self, iatom=None, n_max=3, l_max=6, rcut=None):
        r = self.cutoff if rcut is None else rcut
        if iatom is None:
            return [
                list(spec)
                for spec in _core.soapSpectrumAll(
                    self.cloud, self.neighbor_list, n_max, l_max, r
                )
            ]
        return list(
            _core.soapSpectrum(
                self.cloud, iatom, self.neighbor_list, n_max, l_max, r
            )
        )

    def voronoi_features(self, cutoff=None):
        cut = self.cutoff if cutoff is None else cutoff
        return [list(row) for row in _core.voronoiFeatures(self.cloud, cut)]

    def fit_classifier(self, X, y, labels=None):
        clf = _core.LinearClassifier()
        if labels is not None:
            clf.labels = list(labels)
        clf.fit(X, y)
        self._classifier = clf
        return clf

    def predict_class(self, x):
        if self._classifier is None:
            raise RuntimeError("fit_classifier must be called first")
        return self._classifier.predict(x)


def read(filename, **kwargs):
    """Load a LAMMPS dump. See Frame.from_file."""
    return Frame.from_file(filename, **kwargs)
