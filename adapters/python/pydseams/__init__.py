# Copyright (c) d-SEAMS developers
"""d-SEAMS Python front end with multi-frame TrajectorySequence and time-correlation tools."""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

_THIS_DIR = Path(__file__).resolve().parent


def _load_downstream_pydseams():
    """Locate and execute an installed site-packages pydseams if present."""
    for entry in sys.path:
        try:
            base = Path(entry).resolve() if entry else Path.cwd().resolve()
        except OSError:
            continue
        cand = base / "pydseams"
        init_py = cand / "__init__.py"
        if cand != _THIS_DIR and init_py.is_file():
            if str(cand) not in __path__:
                __path__.append(str(cand))
            spec = importlib.util.spec_from_file_location(
                __name__,
                str(init_py),
                submodule_search_locations=list(__path__),
            )
            if spec is not None and spec.loader is not None:
                mod = importlib.util.module_from_spec(spec)
                mod.__package__ = __name__
                mod.__path__ = __path__
                spec.loader.exec_module(mod)
                return mod
    return None


_downstream = _load_downstream_pydseams()
if _downstream is not None:
    for _k, _v in vars(_downstream).items():
        if not _k.startswith("__"):
            globals()[_k] = _v
    __version__ = getattr(_downstream, "__version__", "2.14.3")
else:
    from .trajectory import CtypesFrame as Frame, _load_frame_at_step

    __version__ = "2.14.3"

    def read(path, frame=1, **kwargs):
        """Load one configuration via libseams."""
        return _load_frame_at_step(Path(path), int(frame), **kwargs)

    def from_arrays(positions, cell, numbers=None, cutoff=None, bonded="cutoff"):
        """Construct a single Frame from (N, 3) positions and a periodic cell."""
        return Frame(positions, cell, types=numbers, cutoff=3.5 if cutoff is None else cutoff)


from .trajectory import (  # noqa: E402
    DiffusionFit,
    HbondLifetimeAccumulator,
    HbondLifetimeFit,
    HbondLifetimeResult,
    MoleculeCentresResult,
    MoleculeDipoles,
    MsdAccumulator,
    MsdResult,
    OrientationOrder,
    ReorientationAccumulator,
    ReorientationFit,
    ReorientationResult,
    TrajectorySequence,
    UnwrappedCentresResult,
    fit_diffusion,
    fit_exponential_tau,
    hbond_lifetime,
    molecule_centres,
    molecule_dipoles,
    orientation_order,
    reorientation_tcf,
    unwrap_centres,
    unwrap_positions,
    water_hbond_edges,
    water_vectors,
)

# Trajectory stays the single-frame Frame alias for all of 2.14.x; a future
# major version may rebind it to TrajectorySequence.
Trajectory = Frame
open_trajectory = TrajectorySequence.open

_EXTRA_EXPORTS = [
    "DiffusionFit",
    "HbondLifetimeAccumulator",
    "HbondLifetimeFit",
    "HbondLifetimeResult",
    "MoleculeCentresResult",
    "MoleculeDipoles",
    "MsdAccumulator",
    "MsdResult",
    "OrientationOrder",
    "ReorientationAccumulator",
    "ReorientationFit",
    "ReorientationResult",
    "Trajectory",
    "TrajectorySequence",
    "UnwrappedCentresResult",
    "fit_diffusion",
    "fit_exponential_tau",
    "hbond_lifetime",
    "molecule_centres",
    "molecule_dipoles",
    "open_trajectory",
    "orientation_order",
    "reorientation_tcf",
    "unwrap_centres",
    "unwrap_positions",
    "water_hbond_edges",
    "water_vectors",
]

_FRAME_EXPORTS = ("Frame", "read", "from_arrays", "__version__")
__all__ = sorted(set(getattr(_downstream, "__all__", _FRAME_EXPORTS)).union(_EXTRA_EXPORTS))
