# Copyright (c) d-SEAMS developers
"""Status codes and the exception raised from ``seams_handle_error``."""

from __future__ import annotations


OK = 0
ERR_NULL = 1
ERR_RANGE = 2
ERR_BOX = 3
ERR_ARG = 4
ERR_STATE = 5
ERR_NOMEM = 6
ERR_INTERNAL = 7


class SeamsError(Exception):
    """A failing ``libseams`` call.

    Parameters
    ----------
    status : int
        ``seams_status`` value. Zero is never stored here.
    message : str
        Text from ``seams_handle_error``. The library owns that string.
        This object keeps a Python copy.

    Examples
    --------
    >>> err = SeamsError(5, "no frame")
    >>> err.status
    5
    >>> "no frame" in str(err)
    True
    """

    def __init__(self, status: int, message: str) -> None:
        self.status = status
        self.message = message
        super().__init__(f"{message} ({status})")
