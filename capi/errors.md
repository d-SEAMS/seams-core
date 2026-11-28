<!-- https://docs.dseams.info/capi/errors.md -->

# Errors {#errors}

Status codes, the thread-local last error, and version strings.

The stable C interface provides 120 exported functions through `include/seams.h`.

```c
// Copyright (c) d-SEAMS developers
/* A missing path stores a thread-local error and a non-OK status. */
#include "seams.h"

#include <stdio.h>
#include <string.h>

int main(void) {
  seams_system *sys = seams_system_create();
  const char *null_message = NULL;
  seams_status missing =
      seams_system_natoms(NULL, NULL);
  null_message = seams_last_error_message();
  if (missing != SEAMS_ERR_NULL || null_message == NULL || null_message[0] == '\0' ||
      strcmp(null_message, seams_last_error()) != 0) {
    fprintf(stderr, "null status %d message %s\n", (int)missing,
            null_message ? null_message : "(null)");
    return 1;
  }
  seams_status status =
      seams_system_read(sys, "input/traj/no-such-frame.lammpstrj", 1, -1);
  const char *err = seams_last_error();
  const char *message = seams_last_error_message();
  const char *name = seams_status_string(status);
  printf("status %d\n", (int)status);
  if (status == SEAMS_OK || err == NULL || err[0] == '\0' || message == NULL ||
      strcmp(err, message) != 0 || name == NULL || name[0] == '\0') {
    fprintf(stderr, "status %d err %s\n", (int)status, err ? err : "(null)");
    return 1;
  }
  seams_system_destroy(sys);
  return 0;
}
```

## seams_abi {#seams-abi}

```c
int seams_abi(void)
```

ABI number compiled into this library.

Returns: SEAMS_ABI, which is 1.

@par Thread safety Pure. @par Example @code if (seams_abi() != SEAMS_ABI) { return 1; } @endcode

## seams_version {#seams-version}

```c
int seams_version(void)
```

Packed release number.

Returns: SEAMS_VERSION_NUMBER, (major << 16) | (minor << 8) | patch.

@par Thread safety Pure. @par Example @code if (seams_version() < ((2 << 16) | (14 << 8))) { return 1; } @endcode

## seams_version_major {#seams-version-major}

```c
int seams_version_major(void)
```

Major version of the library at run time.

Returns: SEAMS_VERSION_MAJOR of the build.

@par Thread safety Pure. @par Example @code int major = seams_version_major(); @endcode

## seams_version_minor {#seams-version-minor}

```c
int seams_version_minor(void)
```

Minor version of the library at run time.

Returns: SEAMS_VERSION_MINOR of the build.

@par Thread safety Pure. @par Example @code int minor = seams_version_minor(); @endcode

## seams_version_patch {#seams-version-patch}

```c
int seams_version_patch(void)
```

Patch version of the library at run time.

Returns: SEAMS_VERSION_PATCH of the build.

@par Thread safety Pure. @par Example @code int patch = seams_version_patch(); @endcode

## seams_version_string {#seams-version-string}

```c
const char *seams_version_string(void)
```

Version text of the built library.

Returns: Pointer to a static string. Never NULL.

@par Ownership The caller does not free it. @par Thread safety Pure. @par Example @code printf(“d-SEAMS %s\\n”, seams_version_string()); @endcode

## seams_status_string {#seams-status-string}

```c
const char *seams_status_string(seams_status status)
```

Short name of a status code.

**`status`**

A ::seams_status value.

Returns: Pointer to a static string. Unknown codes return “unknown status”.

@par Ownership The caller does not free it. @par Thread safety Pure. @par Example @code fprintf(stderr, “%s\\n”, seams_status_string(SEAMS_ERR_BOX)); @endcode

## seams_last_error {#seams-last-error}

```c
const char *seams_last_error(void)
```

Thread-local detail of the last failing call on this thread.

Returns: Pointer to a thread-local string. Never NULL.

ABI 1 spelling of seams_last_error_message(). Both return the same string.

@par Ownership Owned by the library. Valid until the next call on this thread that sets or clears the error. @par Thread safety The storage is thread-local. @par Example @code if (seams_system_set_atoms(sys, xyz, n, box9, 9, NULL, NULL) != SEAMS_OK) { fprintf(stderr, “%s\\n”, seams_last_error()); } @endcode

## seams_handle_error {#seams-handle-error}

```c
const char *seams_handle_error(const void *handle)
```

Detail string stored on a handle.

**`handle`**

A ::seams_system, ::seams_neighbours, ::seams_rings, ::seams_cage_set, ::seams_topo, or ::seams_brine, or NULL.

Returns: The handle’s string, or the thread-local string when handle is NULL. Never NULL.

A call taking the handle non-const clears the string on entry. A call taking it const writes the string only when it fails, so after a successful const call the string still holds the last failure. @par Ownership Owned by the handle. Valid until the next failing call on that handle, the next non-const call on it, or destroy. @par Thread safety Reading a handle requires that no other thread is inside a call on that same handle. @par Example @code if (seams_chill_plus_n(sys, 3.5, 0, SEAMS_GRAPH_CUTOFF, SEAMS_CHILL_PUBLISHED, labels, n, NULL) != SEAMS_OK) { fprintf(stderr, “%s\\n”, seams_handle_error(sys)); } @endcode

## seams_label_name {#seams-label-name}

```c
const char *seams_label_name(int label)
```

Name of a CHILL label.

**`label`**

A ::seams_label value.

Returns: Static string. Unknown codes return “unknown”.

@par Ownership The caller does not free it. @par Thread safety Pure. @par Example @code printf(“%s\\n”, seams_label_name(SEAMS_LABEL_HEXAGONAL)); @endcode

## seams_set_num_threads {#seams-set-num-threads}

```c
seams_status seams_set_num_threads(int n)
```

Set the OpenMP team used by later analysis calls.

**`n`**

Positive width, or 0 to restore the width captured at first use. A negative width is rejected.

Returns: ::SEAMS_OK, or ::SEAMS_ERR_RANGE when n is negative, ::SEAMS_ERR_NOMEM, ::SEAMS_ERR_INTERNAL.

@par Thread safety Process-global. Do not call it concurrently with analysis. The stored width is then read by every handle. @par Example @code seams_set_num_threads(1); @endcode

## seams_get_num_threads {#seams-get-num-threads}

```c
int seams_get_num_threads(void)
```

Team width the next analysis call will request.

Returns: A positive count. After seams_set_num_threads(0) this is the saved default, which is omp_get_max_threads() when OpenMP is linked and 1 otherwise.

@par Thread safety Observes the process-global width. Do not call it concurrently with seams_set_num_threads(). @par Example @code int width = seams_get_num_threads(); @endcode

## seams_last_error_message {#seams-last-error-message}

```c
const char *seams_last_error_message(void)
```

Thread-local detail of the last failing call on this thread.

Returns: Pointer to a thread-local string. Never NULL.

The string is empty after a successful call on this thread. Prefer seams_handle_error() when the failing call had a handle: that copy survives a later call on another handle. seams_last_error() is the ABI 1 spelling and returns the same string.

@par Ownership Owned by the library. Valid until the next call on this thread that sets or clears the error. @par Thread safety The storage is thread-local. @par Example @code seams_status rc = seams_system_natoms(NULL, NULL); fprintf(stderr, “%s: %s\\n”, seams_status_string(rc), seams_last_error_message()); @endcode
