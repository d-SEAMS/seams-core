R
=

roxygen2 comments in ``bindings/r/dseams/R/dseams.R`` are the source of
``man/*.Rd``. Every ``NAMESPACE`` export has an ``\\alias``.
``R CMD check --as-cran`` must not report an undocumented object.

The ``chill_plus`` example reads Ice Ih when ``DSEAMS_ROOT`` is set and
asserts 128 hexagonal waters and 256 six-rings:

.. code-block:: r

   frame <- read_frame(file.path(root, "input", "traj", "ice1h.lammpstrj"),
                       type = -1L)
   ice <- chill_plus(frame, graph = "cutoff")
   stopifnot(sum(ice$name == "hexagonal") == 128L)
   rings <- primitive_rings(frame, graph = "cutoff", max_depth = 6L)
   stopifnot(rings$n[rings$size == 6L] == 256L)
