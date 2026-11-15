# Copyright (c) d-SEAMS developers
"""ctypes bindings for the d-SEAMS C API.

``import seams`` loads no shared library. The first call resolves
``SEAMS_LIB`` or the loader name ``libseams``.

Three Python surfaces read the same engine:

* ``seams`` is this package. It wraps ``include/seams.h`` with NumPy
  arrays, dataclasses, and :class:`seams.SeamsError`.
* ``dseams`` is the in-tree nanobind module. It covers ``chill_plus``,
  ``q4q6``, neighbours, rings, and :class:`dseams.Frame`.
* ``pydseams`` is ``pydseamslib`` on PyPI. It is the trajectory
  :class:`pydseams.Frame`. :func:`from_pydseams` copies one of those
  frames into a :class:`System`.

Examples
--------
>>> from seams import as_box9
>>> float(as_box9((1.0, 2.0, 3.0))[3])
1.0
"""

from __future__ import annotations

from seams.errors import (
    ERR_ARG,
    ERR_BOX,
    ERR_INTERNAL,
    ERR_NOMEM,
    ERR_NULL,
    ERR_RANGE,
    ERR_STATE,
    OK,
    SeamsError,
)
from seams.results import (
    BrineFronts,
    BrinePockets,
    BrineProfile,
    BrineRejection,
    CageCounts,
    ChillResult,
    Contacts,
    Density,
    Domains,
    NeighbourCSR,
    OrderResult,
    Prisms,
    Q4Q6,
    Q6Bar,
    Radial,
    RingLabels,
    RunningCN,
    Signature,
)
from seams.system import (
    CAGE_BOTH,
    CAGE_CUBIC,
    CAGE_HEX,
    CAGE_WATER,
    CUTOFF,
    LABEL_CLATHRATE,
    LABEL_CUBIC,
    LABEL_HEXAGONAL,
    LABEL_INTERCLATHRATE,
    LABEL_INTERFACIAL,
    LABEL_RECUBIC,
    LABEL_REHEX,
    LABEL_UNCLASSIFIED,
    LABEL_WATER,
    MUTUAL,
    OVITO,
    PUBLISHED,
    SEEDED,
    Brine,
    CageSet,
    Fingerprint,
    Neighbours,
    Rings,
    System,
    abi,
    as_box9,
    from_pydseams,
    label_name,
    nframes,
    num_threads,
    set_num_threads,
    status_string,
    version,
    version_tuple,
)

__all__ = [
    "CAGE_BOTH",
    "CAGE_CUBIC",
    "CAGE_HEX",
    "CAGE_WATER",
    "CUTOFF",
    "ERR_ARG",
    "ERR_BOX",
    "ERR_INTERNAL",
    "ERR_NOMEM",
    "ERR_NULL",
    "ERR_RANGE",
    "ERR_STATE",
    "LABEL_CLATHRATE",
    "LABEL_CUBIC",
    "LABEL_HEXAGONAL",
    "LABEL_INTERCLATHRATE",
    "LABEL_INTERFACIAL",
    "LABEL_RECUBIC",
    "LABEL_REHEX",
    "LABEL_UNCLASSIFIED",
    "LABEL_WATER",
    "MUTUAL",
    "OK",
    "OVITO",
    "PUBLISHED",
    "SEEDED",
    "Brine",
    "BrineFronts",
    "BrinePockets",
    "BrineProfile",
    "BrineRejection",
    "CageCounts",
    "CageSet",
    "ChillResult",
    "Contacts",
    "Density",
    "Domains",
    "Fingerprint",
    "NeighbourCSR",
    "Neighbours",
    "OrderResult",
    "Prisms",
    "Q4Q6",
    "Q6Bar",
    "Radial",
    "RingLabels",
    "Rings",
    "RunningCN",
    "SeamsError",
    "Signature",
    "System",
    "abi",
    "as_box9",
    "from_pydseams",
    "label_name",
    "nframes",
    "num_threads",
    "set_num_threads",
    "status_string",
    "version",
    "version_tuple",
]
