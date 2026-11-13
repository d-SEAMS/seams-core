C API
=====


Every function exported by ``src/libseams.map`` (the ``seams_*`` version script) on ``capi/polish-2.14``. The concrete list is the 77 ``SEAMS_API`` functions in ``include/seams.h`` at ``0127f41e``. ``meson test docs_capi_symbols`` fails if one of them is absent from these pages.

.. toctree::
   :maxdepth: 1

   atoms
   neighbours
   order
   rings
   cages
   errors
   further
   merge

