Atoms in
========


Read a frame into an opaque handle.

Generated from the comments on the ``SEAMS_API`` declarations in ``include/seams.h``. This tree has 69 exports. ``meson test docs_capi_symbols`` fails if a declaration is missing from these pages, and ``meson test docs_capi_examples`` compiles the example when this page has one.

.. literalinclude:: ../../../docs/examples/capi_atoms.c
   :language: c

.. _capi-seams_system_create:

seams_system_create
-------------------


.. code-block:: c

   seams_system *seams_system_create(void)

Allocate an empty frame.

Returns: A new handle, or NULL on allocation failure. seams_last_error() then explains the failure. The caller owns the handle and passes it to seams_system_destroy().

Thread safety: safe. The new handle is not shared yet.

.. _capi-seams_system_destroy:

seams_system_destroy
--------------------


.. code-block:: c

   void seams_system_destroy(seams_system *sys)

Free a frame.

``sys``
    Handle from seams_system_create(), or NULL.

Ownership: the handle and the frame it stored are released. Neighbour, ring, cage, and brine handles built from it stay valid: they hold their own copies. Thread safety: the caller must not use sys after this returns, and must not call it concurrently with any use of sys.

.. _capi-seams_system_set_atoms:

seams_system_set_atoms
----------------------


.. code-block:: c

   seams_status seams_system_set_atoms(seams_system *sys, const double *xyz, size_t n, const double *box9, size_t box_len, const int *types, const int *ids)

Store points and a general cell. Lengths are size_t.

``sys``
    Frame handle.

``xyz``
    n times 3 Cartesian coordinates, Angstrom.

``n``
    Atom count. Zero is rejected. A count that does not fit in the frame's integer width is rejected.

``box9``
    Cell numbers. box_len selects how many are read.

``box_len``
    Must be 9.

``types``
    Atom types, or NULL to store type 1 for every atom.

``ids``
    Atom ids, or NULL to store 1..n. Ids must be unique.

Returns: SEAMS_OK, or SEAMS_ERR_NULL, SEAMS_ERR_RANGE, SEAMS_ERR_BOX, SEAMS_ERR_ARG, SEAMS_ERR_NOMEM, SEAMS_ERR_INTERNAL.

Ownership: the library copies the arrays. The caller keeps them. Thread safety: exclusive use of sys.

.. _capi-seams_system_set_frame:

seams_system_set_frame
----------------------


.. code-block:: c

   seams_status seams_system_set_frame(seams_system *sys, const double *xyz, int n, const double *box9, const int *types, const int *ids)

Store points and a general LAMMPS cell.

``sys``
    Frame handle.

``xyz``
    n times 3 coordinates, Angstrom.

``n``
    Atom count. A non-positive count is rejected.

``box9``
    Nine cell numbers: ox, oy, oz, lx, ly, lz, xy, xz, yz.

``types``
    Types, or NULL for type 1.

``ids``
    Ids, or NULL for 1..n.

Returns: Same codes as seams_system_set_atoms().

ABI 1 entry. n is an int. New callers use seams_system_set_atoms(). Ownership: the library copies the arrays. Thread safety: exclusive use of sys.

.. _capi-seams_system_set_ortho:

seams_system_set_ortho
----------------------


.. code-block:: c

   seams_status seams_system_set_ortho(seams_system *sys, const double *xyz, int n, double ox, double oy, double oz, double lx, double ly, double lz, const int *types, const int *ids)

Store points and an orthorhombic cell.

``sys``
    Frame handle.

``xyz``
    n times 3 coordinates, Angstrom.

``n``
    Atom count.

``ox``
    Origin x, Angstrom.

``oy``
    Origin y, Angstrom.

``oz``
    Origin z, Angstrom.

``lx``
    Edge length x, Angstrom.

``ly``
    Edge length y, Angstrom.

``lz``
    Edge length z, Angstrom.

``types``
    Types, or NULL.

``ids``
    Ids, or NULL.

Returns: Same codes as seams_system_set_frame(). The tilt is zero.

Thread safety: exclusive use of sys.

.. _capi-seams_system_read:

seams_system_read
-----------------


.. code-block:: c

   seams_status seams_system_read(seams_system *sys, const char *path, int frame, int type_filter)

Read one frame from a trajectory file.

``sys``
    Frame handle. On failure the previously stored frame is kept.

``path``
    Filesystem path. CON, LAMMPS dump, and the formats the linked readers accept.

``frame``
    Frame index, 1-based, matching the file order.

``type_filter``
    Negative keeps every atom. Zero guesses the oxygen species. A positive value keeps that species number.

Returns: SEAMS_OK, or SEAMS_ERR_NULL, SEAMS_ERR_RANGE, SEAMS_ERR_ARG, SEAMS_ERR_NOMEM, SEAMS_ERR_INTERNAL.

Ownership: the library copies the frame out of the file. Thread safety: exclusive use of sys. Concurrent reads of distinct paths on distinct handles are safe.

.. _capi-seams_path_nframes:

seams_path_nframes
------------------


.. code-block:: c

   seams_status seams_path_nframes(const char *path, size_t *nframes)

Number of frames in a trajectory path.

``path``
    Filesystem path.

``nframes``
    Receives the count. Must be non-NULL.

Returns: SEAMS_OK, or SEAMS_ERR_NULL, SEAMS_ERR_RANGE, SEAMS_ERR_ARG, SEAMS_ERR_INTERNAL. \*nframes is unchanged on failure.

Thread safety: safe for distinct paths.

.. _capi-seams_system_count:

seams_system_count
------------------


.. code-block:: c

   int seams_system_count(const seams_system *sys)

Number of points in the frame.

``sys``
    Frame handle, or NULL.

Returns: The atom count, or 0 when sys is NULL or holds no frame.

Thread safety: exclusive use of sys.

