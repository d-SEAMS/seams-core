Errors
======


Status codes, the thread-local last error, and version strings.

Generated from the comments on the ``SEAMS_API`` declarations in ``include/seams.h``. This tree has 83 exports. ``meson test docs_capi_symbols`` fails if a declaration is missing from these pages, and ``meson test docs_capi_examples`` compiles the example when this page has one.

.. literalinclude:: ../../../docs/examples/capi_errors.c
   :language: c

.. _capi-seams_abi:

seams_abi
---------


.. code-block:: c

   int seams_abi(void)

ABI number compiled into this library.

Returns: SEAMS_ABI, which is 1.

@par Thread safety Pure. @par Example @code if (seams_abi() != SEAMS_ABI) { return 1; } @endcode

.. _capi-seams_version:

seams_version
-------------


.. code-block:: c

   int seams_version(void)

Packed release number.

Returns: SEAMS_VERSION_NUMBER, (major << 16) | (minor << 8) | patch.

@par Thread safety Pure. @par Example @code if (seams_version() < ((2 << 16) | (14 << 8))) { return 1; } @endcode

.. _capi-seams_version_major:

seams_version_major
-------------------


.. code-block:: c

   int seams_version_major(void)

Major version of the library at run time.

Returns: SEAMS_VERSION_MAJOR of the build.

@par Thread safety Pure. @par Example @code int major = seams_version_major(); @endcode

.. _capi-seams_version_minor:

seams_version_minor
-------------------


.. code-block:: c

   int seams_version_minor(void)

Minor version of the library at run time.

Returns: SEAMS_VERSION_MINOR of the build.

@par Thread safety Pure. @par Example @code int minor = seams_version_minor(); @endcode

.. _capi-seams_version_patch:

seams_version_patch
-------------------


.. code-block:: c

   int seams_version_patch(void)

Patch version of the library at run time.

Returns: SEAMS_VERSION_PATCH of the build.

@par Thread safety Pure. @par Example @code int patch = seams_version_patch(); @endcode

.. _capi-seams_version_string:

seams_version_string
--------------------


.. code-block:: c

   const char *seams_version_string(void)

Version text, for example "2.14.0".

Returns: Pointer to a static string. Never NULL.

@par Ownership The caller does not free it. @par Thread safety Pure. @par Example @code printf("d-SEAMS %s\\n", seams_version_string()); @endcode

.. _capi-seams_status_string:

seams_status_string
-------------------


.. code-block:: c

   const char *seams_status_string(seams_status status)

Short name of a status code.

``status``
    A ::seams_status value.

Returns: Pointer to a static string. Unknown codes return "unknown status".

@par Ownership The caller does not free it. @par Thread safety Pure. @par Example @code fprintf(stderr, "%s\\n", seams_status_string(SEAMS_ERR_BOX)); @endcode

.. _capi-seams_last_error:

seams_last_error
----------------


.. code-block:: c

   const char *seams_last_error(void)

Thread-local detail of the last failing call on this thread.

Returns: Pointer to a thread-local string. Never NULL.

ABI 1 spelling of seams_last_error_message(). Both return the same string.

@par Ownership Owned by the library. Valid until the next call on this thread that sets or clears the error. @par Thread safety The storage is thread-local. @par Example @code if (seams_system_set_atoms(sys, xyz, n, box9, 9, NULL, NULL) != SEAMS_OK) { fprintf(stderr, "%s\\n", seams_last_error()); } @endcode

.. _capi-seams_handle_error:

seams_handle_error
------------------


.. code-block:: c

   const char *seams_handle_error(const void *handle)

Detail string stored on a handle.

``handle``
    A ::seams_system, ::seams_neighbours, ::seams_rings, ::seams_cage_set, ::seams_topo, or ::seams_brine, or NULL.

Returns: The handle's string, or the thread-local string when handle is NULL. Never NULL.

A call taking the handle non-const clears the string on entry. A call taking it const writes the string only when it fails, so after a successful const call the string still holds the last failure. @par Ownership Owned by the handle. Valid until the next failing call on that handle, the next non-const call on it, or destroy. @par Thread safety Reading a handle requires that no other thread is inside a call on that same handle. @par Example @code if (seams_chill_plus_n(sys, 3.5, 0, SEAMS_GRAPH_CUTOFF, SEAMS_CHILL_PUBLISHED, labels, n, NULL) != SEAMS_OK) { fprintf(stderr, "%s\\n", seams_handle_error(sys)); } @endcode

.. _capi-seams_label_name:

seams_label_name
----------------


.. code-block:: c

   const char *seams_label_name(int label)

Name of a CHILL label.

``label``
    A ::seams_label value.

Returns: Static string. Unknown codes return "unknown".

@par Ownership The caller does not free it. @par Thread safety Pure. @par Example @code printf("%s\\n", seams_label_name(SEAMS_LABEL_HEXAGONAL)); @endcode

.. _capi-seams_set_num_threads:

seams_set_num_threads
---------------------


.. code-block:: c

   seams_status seams_set_num_threads(int n)

Set the OpenMP team used by later analysis calls.

``n``
    Positive width, or 0 to restore the width captured at first use. A negative width is rejected.

Returns: ::SEAMS_OK, or ::SEAMS_ERR_RANGE when n is negative, ::SEAMS_ERR_NOMEM, ::SEAMS_ERR_INTERNAL.

@par Thread safety Process-global. Do not call it concurrently with analysis. The stored width is then read by every handle. @par Example @code seams_set_num_threads(1); @endcode

.. _capi-seams_get_num_threads:

seams_get_num_threads
---------------------


.. code-block:: c

   int seams_get_num_threads(void)

Team width the next analysis call will request.

Returns: A positive count. After seams_set_num_threads(0) this is the saved default, which is omp_get_max_threads() when OpenMP is linked and 1 otherwise.

@par Thread safety Observes the process-global width. Do not call it concurrently with seams_set_num_threads(). @par Example @code int width = seams_get_num_threads(); @endcode

.. _capi-seams_last_error_message:

seams_last_error_message
------------------------


.. code-block:: c

   const char *seams_last_error_message(void)

Thread-local detail of the last failing call on this thread.

Returns: Pointer to a thread-local string. Never NULL.

The string is empty after a successful call on this thread. Prefer seams_handle_error() when the failing call had a handle: that copy survives a later call on another handle. seams_last_error() is the ABI 1 spelling and returns the same string.

@par Ownership Owned by the library. Valid until the next call on this thread that sets or clears the error. @par Thread safety The storage is thread-local. @par Example @code seams_status rc = seams_system_natoms(NULL, NULL); fprintf(stderr, "%s: %s\\n", seams_status_string(rc), seams_last_error_message()); @endcode

