from . import _core
from pathlib import Path
import tempfile

__version__ = "2.0.0"


class Trajectory:
    """High-level interface for d-SEAMS ice structure analysis."""

    def __init__(self, filename, frame=1, atom_type=2, cutoff=3.5):
        """Load a LAMMPS trajectory frame.

        Parameters
        ----------
        filename : str or Path
            Path to LAMMPS trajectory file.
        frame : int
            Frame number to read (1-indexed).
        atom_type : int
            Type ID of atoms to analyze (default 2 for oxygen).
        cutoff : float
            Neighbor list cutoff distance in Angstroms.
        """
        self.filename = str(Path(filename).resolve())
        self.frame = frame
        self.atom_type = atom_type
        self.cutoff = cutoff

        self.cloud = _core.readLammpsTrjreduced(
            filename=self.filename,
            targetFrame=frame,
            typeI=atom_type,
            isSlice=False,
            coordLow=[0, 0, 0],
            coordHigh=[0, 0, 0],
        )
        self._nlist = None
        self._hbonds = None
        self._rings = None
        self._classifier = None
        self._ring_updater = _core.RingUpdater(6)

    @property
    def n_atoms(self):
        """Number of atoms in the frame."""
        return self.cloud.nop

    @property
    def box(self):
        """Simulation box dimensions [x, y, z]."""
        return list(self.cloud.box)

    @property
    def neighbor_list(self):
        """Neighbor list (computed lazily)."""
        if self._nlist is None:
            self._nlist = _core.neighListO(
                rcutoff=self.cutoff, yCloud=self.cloud, typeI=self.atom_type
            )
        return self._nlist

    @property
    def hbonds(self):
        """Hydrogen bond network (computed lazily)."""
        if self._hbonds is None:
            self._hbonds = _core.populateHbonds(
                filename=self.filename,
                yCloud=self.cloud,
                nList=self.neighbor_list,
                targetFrame=self.frame,
                Htype=1,
            )
        return self._hbonds

    def load_frame(self, frame):
        """Read another frame, keeping the incremental ring updater."""
        self.frame = frame
        self.cloud = _core.readLammpsTrjreduced(
            filename=self.filename,
            targetFrame=frame,
            typeI=self.atom_type,
            isSlice=False,
            coordLow=[0, 0, 0],
            coordHigh=[0, 0, 0],
        )
        self._nlist = None
        self._hbonds = None
        self._rings = None

    @property
    def rings(self):
        """Primitive ring network (max depth 6).

        Consecutive load_frame calls re-enumerate only sources inside the
        locality bound of a changed bond. The first access is a full pass.
        """
        if self._rings is None:
            hbonds_idx = _core.neighbourListByIndex(
                yCloud=self.cloud, nList=self.hbonds
            )
            self._rings = self._ring_updater.update(hbonds_idx)
        return self._rings

    @property
    def rings_recomputed_sources(self):
        """Sources re-enumerated by the last rings update."""
        return self._ring_updater.lastRecomputedSources()

    def classify_chill_plus(self):
        """Run CHILL+ classification on all atoms.

        Returns
        -------
        dict
            Counts of each ice type: cubic, hexagonal, interfacial,
            clathrate, interClathrate, water, unclassified.
        """
        # Step 1: compute bond correlations
        _core.getCorrelPlus(
            yCloud=self.cloud, nList=self.neighbor_list, isSlice=False)
        # Step 2: classify ice types based on correlations
        outdir = tempfile.mkdtemp(prefix="dseams_")
        _core.getIceTypePlus(
            yCloud=self.cloud, nList=self.neighbor_list,
            path=outdir + "/", firstFrame=self.frame,
            isSlice=False, outputFileName="chillPlus.txt"
        )

        counts = {}
        for pt in self.cloud.pts:
            name = pt.iceType.name
            counts[name] = counts.get(name, 0) + 1
        return counts

    def classify_chill(self):
        """Run CHILL classification on all atoms.

        Returns
        -------
        dict
            Counts of each ice type.
        """
        # Step 1: compute bond correlations
        _core.getCorrel(
            yCloud=self.cloud, nList=self.neighbor_list, isSlice=False)
        # Step 2: classify ice types based on correlations
        outdir = tempfile.mkdtemp(prefix="dseams_")
        _core.getIceType(
            yCloud=self.cloud, nList=self.neighbor_list,
            path=outdir + "/", firstFrame=self.frame,
            isSlice=False, outputFileName="chill.txt"
        )

        counts = {}
        for pt in self.cloud.pts:
            name = pt.iceType.name
            counts[name] = counts.get(name, 0) + 1
        return counts

    def find_prisms(self, output_dir="output/", max_depth=6, shape_matching=False):
        """Run prism analysis for quasi-1D ice.

        Parameters
        ----------
        output_dir : str
            Directory for output files.
        max_depth : int
            Maximum ring size to search.
        shape_matching : bool
            Whether to do shape matching.
        """
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

    def steinhardt(self, order_l=6):
        """Compute the local and neighbour-averaged Steinhardt parameters.

        Parameters
        ----------
        order_l : int
            Degree l of the Steinhardt parameter (3, 4 or 6).

        Returns
        -------
        dict
            "ql" and "ql_bar", each a list of per-atom floats.
        """
        result = _core.steinhardtQl(
            yCloud=self.cloud, nList=self.neighbor_list, orderL=order_l
        )
        return {"ql": list(result.ql), "ql_bar": list(result.qlBar)}

    def steinhardt_voronoi(self, order_l=6, cutoff=None):
        """Voronoi facet-area weighted Steinhardt parameters (Mickel)."""
        result = _core.steinhardtQlVoronoi(
            yCloud=self.cloud,
            candidateCutoff=cutoff if cutoff is not None else self.cutoff,
            orderL=order_l,
        )
        return {"ql": list(result.ql), "ql_bar": list(result.qlBar)}

    def classify_templates(self, k_neigh=12):
        """Overlay each neighbour shell onto FCC, HCP, BCC and SC templates.

        Returns
        -------
        list of dict
            One entry per atom with keys "name", "rmsd" and "kind".
            "kind" is the template name string (fcc, hcp, bcc, sc, other).
        """
        hits = _core.classifyTemplates(self.cloud, self.neighbor_list, k_neigh)
        rows = []
        for h in hits:
            name = h.name
            kind = name or getattr(h.kind, "name", str(h.kind))
            rows.append({"name": name, "rmsd": h.rmsd, "kind": kind})
        return rows

    def soap(self, iatom=None, n_max=3, l_max=6, rcut=None):
        """SOAP power spectrum of one particle, or of every atom.

        Parameters
        ----------
        iatom : int or None
            Atom index. None returns one spectrum per atom.
        n_max : int
            Number of radial Gaussians.
        l_max : int
            Highest spherical-harmonic degree.
        rcut : float or None
            SOAP cutoff. Defaults to the Trajectory cutoff.

        Returns
        -------
        list
            A single spectrum when iatom is an int, or a list of spectra
            when iatom is None.
        """
        r = self.cutoff if rcut is None else rcut
        if iatom is None:
            soap_all = getattr(_core, "soapSpectrumAll", None)
            if soap_all is not None:
                return [
                    list(spec)
                    for spec in soap_all(
                        self.cloud, self.neighbor_list, n_max, l_max, r
                    )
                ]
            return [
                list(
                    _core.soapSpectrum(
                        self.cloud, i, self.neighbor_list, n_max, l_max, r
                    )
                )
                for i in range(self.n_atoms)
            ]
        return list(
            _core.soapSpectrum(
                self.cloud, iatom, self.neighbor_list, n_max, l_max, r
            )
        )

    def voronoi_features(self, cutoff=None):
        """Per-atom Voronoi-weighted [q4, q6, q8].

        Prefers a single cloud-wide _core.voronoiFeatures call when the
        extension exposes it. Otherwise packs three steinhardtQlVoronoi
        results. Does not loop _core.voronoiFeature.
        """
        cut = self.cutoff if cutoff is None else cutoff
        all_fn = getattr(_core, "voronoiFeatures", None)
        if all_fn is not None:
            return [list(row) for row in all_fn(self.cloud, cut)]
        q4 = _core.steinhardtQlVoronoi(
            yCloud=self.cloud, candidateCutoff=cut, orderL=4
        )
        q6 = _core.steinhardtQlVoronoi(
            yCloud=self.cloud, candidateCutoff=cut, orderL=6
        )
        q8 = _core.steinhardtQlVoronoi(
            yCloud=self.cloud, candidateCutoff=cut, orderL=8
        )
        return [
            [q4.ql[i], q6.ql[i], q8.ql[i]] for i in range(self.n_atoms)
        ]

    def fit_classifier(self, X, y, labels=None):
        """Fit a one-versus-rest ridge classifier stored on this trajectory.

        Parameters
        ----------
        X : list of list of float
            Feature matrix, one row per sample.
        y : list of int
            Class indices in ``0 .. n_classes-1``.
        labels : list of str or None
            Optional class names stored on the classifier.

        Returns
        -------
        LinearClassifier
            The fitted _core.LinearClassifier kept on this object.
        """
        clf = _core.LinearClassifier()
        if labels is not None:
            clf.labels = list(labels)
        clf.fit(X, y)
        self._classifier = clf
        return clf

    def predict_class(self, x):
        """Predict the class index of one feature row.

        Parameters
        ----------
        x : list of float
            A single feature row with the same width used in fit_classifier.
        """
        if self._classifier is None:
            raise RuntimeError("fit_classifier must be called first")
        return self._classifier.predict(x)


# Re-export low-level API for advanced users
__all__ = ["Trajectory", "_core", "__version__"]
