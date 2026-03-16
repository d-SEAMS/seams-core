from pydseamslib import _core
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

    @property
    def rings(self):
        """Primitive ring network (computed lazily, max depth 6)."""
        if self._rings is None:
            hbonds_idx = _core.neighbourListByIndex(
                yCloud=self.cloud, nList=self.hbonds
            )
            self._rings = _core.ringNetwork(nList=hbonds_idx, maxDepth=6)
        return self._rings

    def classify_chill_plus(self):
        """Run CHILL+ classification on all atoms.

        Returns
        -------
        dict
            Counts of each ice type: cubic, hexagonal, interfacial,
            clathrate, interClathrate, water, unclassified.
        """
        # Step 1: compute bond correlations
        self.cloud = _core.getCorrelPlus(
            yCloud=self.cloud, nList=self.neighbor_list, isSlice=False)
        # Step 2: classify ice types based on correlations
        outdir = tempfile.mkdtemp(prefix="dseams_")
        self.cloud = _core.getIceTypePlus(
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
        self.cloud = _core.getCorrel(
            yCloud=self.cloud, nList=self.neighbor_list, isSlice=False)
        # Step 2: classify ice types based on correlations
        outdir = tempfile.mkdtemp(prefix="dseams_")
        self.cloud = _core.getIceType(
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


# Re-export low-level API for advanced users
__all__ = ["Trajectory", "_core", "__version__"]
