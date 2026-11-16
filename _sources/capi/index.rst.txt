C API
=====


Every ``SEAMS_API`` function in ``include/seams.h``. The list is generated from that header, so a function that is not declared there is not documented here. ``meson test docs_capi_symbols`` reads the header again and fails if one of those declarations is absent from these pages.

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

