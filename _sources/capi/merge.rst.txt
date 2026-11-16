Where the API lives
===================


``include/seams.h`` in this tree is the C API. ``docs/render_capi.py`` reads that header and writes these pages and ``docs/capi/symbols.txt``.

Regenerate them with:

.. code-block:: bash

   python3 docs/render_capi.py

``meson test docs_capi_symbols`` fails if the header and the pages disagree. ``meson test docs_capi_examples`` compiles and runs the six examples. Those examples call only functions the header declares.
