Errors
======


Status codes, the thread-local last error, and version strings.

Generated from the comments on the ``SEAMS_API`` declarations in ``include/seams.h``. This tree has 69 exports. ``meson test docs_capi_symbols`` fails if a declaration is missing from these pages, and ``meson test docs_capi_examples`` compiles the example when this page has one.

.. literalinclude:: ../../../docs/examples/capi_errors.c
   :language: c

.. _capi-seams_abi:

seams_abi
---------


.. code-block:: c

   int seams_abi(void)

ABI number compiled into this library.

Returns: SEAMS_ABI, which is 1.

Thread safety: pure.

.. _capi-seams_version:

seams_version
-------------


.. code-block:: c

   int seams_version(void)

Packed release number.

Returns: SEAMS_VERSION_NUMBER, (major << 16) | (minor << 8) | patch.

Thread safety: pure.

.. _capi-seams_version_major:

seams_version_major
-------------------


.. code-block:: c

   int seams_version_major(void)

@brief Major version. Thread safety: pure.

.. _capi-seams_version_minor:

seams_version_minor
-------------------


.. code-block:: c

   int seams_version_minor(void)

@brief Minor version. Thread safety: pure.

.. _capi-seams_version_patch:

seams_version_patch
-------------------


.. code-block:: c

   int seams_version_patch(void)

@brief Patch version. Thread safety: pure.

.. _capi-seams_version_string:

seams_version_string
--------------------


.. code-block:: c

   const char *seams_version_string(void)

Version text, for example "2.14.0".

Returns: Pointer to a static string. The caller does not free it.

Thread safety: pure.

.. _capi-seams_status_string:

seams_status_string
-------------------


.. code-block:: c

   const char *seams_status_string(seams_status status)

Short name of a status code.

``status``
    A seams_status value.

Returns: Pointer to a static string. Unknown codes return "unknown status".

Thread safety: pure.

.. _capi-seams_last_error:

seams_last_error
----------------


.. code-block:: c

   const char *seams_last_error(void)

Thread-local detail of the last failing call on this thread.

Returns: Pointer to a thread-local string. Never NULL. The caller does not free it.

The string is empty after a successful call on this thread. It is owned by the library and remains valid until the next call on this thread that sets or clears the error. Prefer seams_handle_error() when the failing call had a handle: that copy survives a later call on another handle.

Thread safety: the storage is thread-local.

.. _capi-seams_handle_error:

seams_handle_error
------------------


.. code-block:: c

   const char *seams_handle_error(const void *handle)

Detail string stored on a handle.

``handle``
    A seams_system, seams_neighbours, seams_rings, seams_cage_set, seams_topo, or seams_brine, or NULL.

Returns: The handle's string, or the thread-local string when handle is NULL. Never NULL. The caller does not free it. The handle string stays valid until the next call on that handle or until destroy.

Thread safety: reading a handle requires that no other thread is inside a call on that same handle.

.. _capi-seams_label_name:

seams_label_name
----------------


.. code-block:: c

   const char *seams_label_name(int label)

Name of a CHILL label.

``label``
    A seams_label value.

Returns: Static string. Unknown codes return "unknown".

Thread safety: pure.

.. _capi-seams_set_num_threads:

seams_set_num_threads
---------------------


.. code-block:: c

   seams_status seams_set_num_threads(int n)

Set the OpenMP team used by later analysis calls.

``n``
    Positive width, or 0 to restore the width captured at first use. A negative width is rejected.

Returns: SEAMS_OK, or SEAMS_ERR_RANGE when n is negative.

Thread safety: process-global. Do not call it concurrently with analysis. The stored width is then read by every handle.

.. _capi-seams_get_num_threads:

seams_get_num_threads
---------------------


.. code-block:: c

   int seams_get_num_threads(void)

Team width the next analysis call will request.

Returns: A positive count. After seams_set_num_threads(0) this is the saved default, which is omp_get_max_threads() when OpenMP is linked and 1 otherwise.

Thread safety: observes the process-global width. Do not call it concurrently with seams_set_num_threads().

