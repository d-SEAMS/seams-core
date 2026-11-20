# Copyright (c) d-SEAMS developers
"""d-SEAMS Python front end with multi-frame Trajectory and time-correlation tools."""

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
                "_pydseams_downstream",
                str(init_py),
                submodule_search_locations=[str(cand)],
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
    __version__ = getattr(_downstream, "__version__", "2.14.2")
else:
    from .trajectory import CtypesFrame as Frame, _load_frame_at_step

    __version__ = "2.14.2"

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
    Trajectory,
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

open_trajectory = Trajectory.open

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

if "__all__" in globals():
    __all__ = sorted(set(globals()["__all__"]).union(_EXTRA_EXPORTS))
else:
    __all__ = sorted(set(_EXTRA_EXPORTS) | {"Frame", "read", "from_arrays", "__version__"})
