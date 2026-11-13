Atoms in
========


Read a frame into an opaque handle.

Generated from the Doxygen XML of ``include/seams.h`` at ``capi/polish-2.14`` commit ``0127f41efc3f44cefdcf532bdea3bba36b05176d``. The C sources of that branch are not part of this docs tree.

``meson test docs_capi_examples`` compiles and runs this example.

.. literalinclude:: ../../../docs/examples/capi_atoms.c
   :language: c

.. _capi-seams_system_create:

seams_system_create
-------------------


.. code-block:: c

   seams_system * seams_system_create(void)

Allocate an empty frame.

A new handle, or NULL on allocation failure. seams_last_error_message() then explains the failure. Ownership The caller owns the handle and passes it to seams_system_destroy(). Thread safety Safe. The new handle is not shared yet. Example seams_system\*sys=seams_system_create(); if(sys==NULL){returnSEAMS_ERR_NOMEM;} seams_system_destroy(sys);

.. _capi-seams_system_destroy:

seams_system_destroy
--------------------


.. code-block:: c

   void seams_system_destroy(seams_system *sys)

Free a frame.

sys Handle from seams_system_create(), or NULL. NULL is ignored. Ownership The handle and the frame it stored are released. Neighbour, ring, cage, and brine handles built from it stay valid: they hold their own copies. Thread safety The caller must not use sys after this returns, and must not call it concurrently with any use of sys. Example seams_system_destroy(sys); sys=NULL;

.. _capi-seams_system_set_atoms:

seams_system_set_atoms
----------------------


.. code-block:: c

   seams_status seams_system_set_atoms(seams_system *sys, const double *xyz, size_t n, const double *box9, size_t box_len, const int *types, const int *ids)

Store points and a general cell. Lengths are size_t.

On failure the previously stored frame, if any, is kept. sys Frame handle. xyz n times 3 Cartesian coordinates, Angstrom. Every value must be finite. n Atom count. Zero is rejected. A count that does not fit in the frame's integer width is rejected. box9 Cell numbers: ox, oy, oz, lx, ly, lz, xy, xz, yz. box_len Must be 9. types Atom types, or NULL to store type 1 for every atom. ids Atom ids, or NULL to store 1..n. Ids must be unique. SEAMS_OK; SEAMS_ERR_NULL for a NULL sys, xyz, or box9; SEAMS_ERR_RANGE for n == 0 or box_len != 9; SEAMS_ERR_BOX for a non-finite cell or a non-positive edge; SEAMS_ERR_ARG for a non-finite coordinate or repeated ids; SEAMS_ERR_NOMEM; SEAMS_ERR_INTERNAL. Ownership The library copies the arrays. The caller keeps them. Thread safety Exclusive use of sys. Example constdoublexyz[6]={0,0,0,1.5,0,0}; constdoublebox9[9]={0,0,0,10,10,10,0,0,0}; seams_statusrc=seams_system_set_atoms(sys,xyz,2,box9,9,NULL,NULL);

.. _capi-seams_system_set_frame:

seams_system_set_frame
----------------------


.. code-block:: c

   seams_status seams_system_set_frame(seams_system *sys, const double *xyz, int n, const double *box9, const int *types, const int *ids)

Store points and a general LAMMPS cell, with an int count.

ABI 1 entry. New callers use seams_system_set_atoms(). sys Frame handle. xyz n times 3 coordinates, Angstrom. n Atom count. A non-positive count is rejected. box9 Nine cell numbers: ox, oy, oz, lx, ly, lz, xy, xz, yz. types Types, or NULL for type 1. ids Ids, or NULL for 1..n. Same codes as seams_system_set_atoms(). Ownership The library copies the arrays. Thread safety Exclusive use of sys. Example seams_statusrc=seams_system_set_frame(sys,xyz,2,box9,NULL,NULL);

.. _capi-seams_system_set_ortho:

seams_system_set_ortho
----------------------


.. code-block:: c

   seams_status seams_system_set_ortho(seams_system *sys, const double *xyz, int n, double ox, double oy, double oz, double lx, double ly, double lz, const int *types, const int *ids)

Store points and an orthorhombic cell.

sys Frame handle. xyz n times 3 coordinates, Angstrom. n Atom count. ox Origin x, Angstrom. oy Origin y, Angstrom. oz Origin z, Angstrom. lx Edge length x, Angstrom. ly Edge length y, Angstrom. lz Edge length z, Angstrom. types Types, or NULL. ids Ids, or NULL. Same codes as seams_system_set_frame(). The tilt is zero. Ownership The library copies the arrays. Thread safety Exclusive use of sys. Example seams_statusrc=seams_system_set_ortho(sys,xyz,2,0,0,0,10,10,10,NULL,NULL);

.. _capi-seams_system_read:

seams_system_read
-----------------


.. code-block:: c

   seams_status seams_system_read(seams_system *sys, const char *path, int frame, int type_filter)

Read one frame from a trajectory file.

sys Frame handle. On failure the previously stored frame is kept. path Filesystem path. CON, LAMMPS dump, and the formats the linked readers accept. frame Frame index, 1-based, matching the file order. type_filter Negative keeps every atom. Zero guesses the oxygen species. A positive value keeps that species number. SEAMS_OK, or SEAMS_ERR_NULL, SEAMS_ERR_RANGE (frame index or empty frame), SEAMS_ERR_ARG (unreadable file), SEAMS_ERR_NOMEM, SEAMS_ERR_INTERNAL. Ownership The library copies the frame out of the file. Thread safety Exclusive use of sys. Concurrent reads of distinct paths on distinct handles are safe. Example seams_statusrc=seams_system_read(sys,"ice.lammpstrj",1,-1);

.. _capi-seams_system_count:

seams_system_count
------------------


.. code-block:: c

   int seams_system_count(const seams_system *sys)

Number of points in the frame, as an int.

ABI 1 query. It cannot report an error; seams_system_natoms() can. sys Frame handle, or NULL. The atom count, or 0 when sys is NULL or holds no frame. Thread safety Exclusive use of sys. Example intn=seams_system_count(sys);

.. _capi-seams_system_natoms:

seams_system_natoms
-------------------


.. code-block:: c

   seams_status seams_system_natoms(const seams_system *sys, size_t *n)

Number of points in the frame.

sys Frame handle. n Receives the atom count. SEAMS_OK, SEAMS_ERR_NULL when sys or n is NULL, or SEAMS_ERR_STATE when sys holds no frame. \*n is unchanged on failure. Thread safety Exclusive use of sys. Example size_tn=0; if(seams_system_natoms(sys,&n)==SEAMS_OK){labels=malloc(n\*sizeof\*labels);}

.. _capi-seams_path_nframes:

seams_path_nframes
------------------


.. code-block:: c

   seams_status seams_path_nframes(const char *path, size_t *nframes)

Number of frames in a trajectory path.

path Filesystem path. nframes Receives the count. Must be non-NULL. SEAMS_OK, or SEAMS_ERR_NULL, SEAMS_ERR_RANGE, SEAMS_ERR_ARG, SEAMS_ERR_INTERNAL. \*nframes is unchanged on failure. Thread safety Safe for distinct paths. Example size_tnframes=0; seams_statusrc=seams_path_nframes("ice.lammpstrj",&nframes);

