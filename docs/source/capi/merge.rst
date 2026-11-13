===================
Where the API lives
===================

This section is generated from the Doxygen XML of ``include/seams.h`` on
``capi/polish-2.14`` at ``0127f41efc3f44cefdcf532bdea3bba36b05176d``.
That branch's C sources are not merged into ``docs/release-2.14``.

The version script on that commit, ``src/libseams.map``, exports
``seams_*``. Expanded against the header, that is 77 ``SEAMS_API``
functions. ``docs/capi/symbols.txt`` is that expansion.
``meson test docs_capi_symbols`` fails if any of those names is missing
from these pages.

Eight of the 77 declarations are not in the ``include/seams.h`` on this
docs tip. A program that calls them needs the polished header, and a
library built from that branch, before it will compile:

* ``seams_chill_plus_n``
* ``seams_last_error_message``
* ``seams_neighbours_nrows``
* ``seams_q4q6_n``
* ``seams_rings_histogram_n``
* ``seams_rings_labels_n``
* ``seams_rings_nrings``
* ``seams_system_natoms``

``src/libseams.map`` itself is also absent from this docs tip. The
header that declares those eight functions is ``include/seams.h`` on
``capi/polish-2.14``. ``include/seams_version.h.in`` is included by
that header and is already on this branch.

Merge order:

1. Merge ``capi/polish-2.14`` (at least ``include/seams.h``, the
   definitions of the eight functions, and ``src/libseams.map``) onto
   the library branch.
2. This book already documents that header. Regenerating it is
   ``doxygen`` on that ``include/seams.h``, then
   ``python3 docs/render_capi.py`` on the XML directory.
3. The six examples under ``docs/examples/capi_*.c`` call only
   functions that exist on this docs tip, so
   ``meson test docs_capi_examples`` links the library built here.
   Examples of the eight new entry points wait until step 1.
