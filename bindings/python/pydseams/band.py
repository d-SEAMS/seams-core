# Copyright (c) d-SEAMS developers
"""Primitive-ring census. The count is seams_frame_rings."""


def ring_census(frame, max_size=8):
    """Count primitive rings by size, up to ``max_size``."""
    return frame.ring_census(max_size)
