# Copyright (c) d-SEAMS developers
"""Environment helper for spawning unsanitized host tools from tests.

The sanitizer suite preloads libasan so unsanitized hosts such as Python can
dlopen sanitized libseams. External build and VCS tools (git, ps, meson,
cmake, ninja) are not built with the sanitizer; inheriting the preload makes
their own exit-time allocations or probe binaries abort under LeakSanitizer
or AddressSanitizer. Strip only the preload variables and keep ASAN_OPTIONS,
LSAN_OPTIONS, UBSAN_OPTIONS and everything else intact so sanitized code
under test stays under the sanitizers.
"""

import os
from collections.abc import Mapping

PRELOAD_VARS = ("LD_PRELOAD", "DYLD_INSERT_LIBRARIES", "ASAN_TEST_PRELOAD")


def without_asan_preload(
    base: Mapping[str, str] | None = None,
    extra: Mapping[str, str] | None = None,
) -> dict[str, str]:
    """Return a copy of the environment with sanitizer preload vars removed."""
    env = dict(os.environ if base is None else base)
    for key in PRELOAD_VARS:
        env.pop(key, None)
    if extra:
        env.update(extra)
    return env
