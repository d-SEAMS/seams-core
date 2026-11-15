#ifndef SEAMS_H_
#define SEAMS_H_

/**
 * @file seams.h
 * @brief d-SEAMS stable C API, ABI 1.
 *
 * Opaque handles. A frame is n Cartesian points, in Angstrom, and a
 * LAMMPS cell of nine doubles: ox, oy, oz, lx, ly, lz, xy, xz, yz.
 * Lengths are the edge lengths. A zero tilt is orthorhombic.
 *
 * @par Status and errors
 * Every call that can fail returns a ::seams_status. Zero is success.
 * A failing call stores a detail string for seams_last_error_message()
 * and, when the call had a handle, for seams_handle_error(). The
 * library never aborts or exits the process, and no C++ exception
 * crosses this boundary: every one is caught and mapped to a status.
 * The only functions that do not return a status are the pure
 * metadata queries (version, names), the destroy functions, and the
 * legacy int count queries kept for ABI 1 callers.
 *
 * @par Handles and ownership
 * Every handle comes from a create, build, find, or analyse call and is
 * released by the matching destroy function. Each destroy function
 * accepts NULL and ignores it. A second destroy of the same pointer is
 * detected, records "double destroy", and frees nothing. A handle built
 * from a frame holds its own copy, so the frame may be destroyed first.
 * The library never returns a buffer the caller must free: outputs are
 * written into caller buffers, and the strings it returns are static or
 * owned by the library.
 *
 * @par Buffers
 * Lengths are size_t. A NULL output buffer with a non-NULL length
 * pointer queries the required count and returns ::SEAMS_OK. A short
 * non-NULL buffer returns ::SEAMS_ERR_RANGE and writes nothing into that
 * buffer. The length pointer still receives the required count.
 *
 * @par Const correctness
 * Input arrays are const. A function that takes a const handle is a
 * read-only accessor: on success it writes nothing to the handle, and on
 * failure it writes only the handle's error slot, under a lock. The
 * frame accessors seams_system_count(), seams_system_natoms(),
 * seams_system_box(), seams_system_positions(), seams_system_types()
 * and seams_system_ids(), and every query on a neighbour, ring, cage,
 * brine, or fingerprint handle, are of this kind. A function that takes
 * a non-const ::seams_system says why in its "Handle access" paragraph.
 * There are three reasons. A load replaces the frame. seams_chill(),
 * seams_chill_plus() and seams_chill_plus_n() write their labels into
 * the frame. Every other analysis clears the handle's error slot on
 * entry, so seams_handle_error() reports that call and no earlier one;
 * that is the only handle state those calls write. None of them changes
 * the stored coordinates, cell, types, or ids.
 *
 * @par Thread safety
 * Distinct handles may be used from distinct threads at the same time.
 * Const-handle calls on one handle may run on several threads at once,
 * provided no call taking that handle non-const runs at the same time.
 * Two non-const calls on one handle must not overlap. When threads share
 * a handle, read failures with seams_last_error_message(), which is
 * thread-local; the handle's slot holds whichever failure landed last.
 * seams_set_num_threads() is process-global.
 *
 * @par ABI policy
 * SEAMS_ABI stays 1 while every shipped function keeps its signature and
 * every shipped enumerator keeps its numeric value. New functions and
 * new enumerators are appended. No public struct carries fields. A
 * later struct, if one is added, starts with a uint32_t size set by the
 * caller to the size of the struct the caller compiled, and the library
 * reads only the prefix that size covers. Only symbols marked
 * ::SEAMS_API are exported; the library is built with hidden visibility.
 *
 * @par Example
 * @code
 * #include <seams.h>
 * #include <stdio.h>
 *
 * int classify(const double *xyz, size_t n, const double box9[9], int *labels) {
 *   seams_system *sys = seams_system_create();
 *   seams_status rc = sys == NULL ? SEAMS_ERR_NOMEM
 *                                 : seams_system_set_atoms(sys, xyz, n, box9, 9, NULL, NULL);
 *   if (rc == SEAMS_OK) {
 *     rc = seams_chill_plus_n(sys, 3.5, 0, SEAMS_GRAPH_CUTOFF, SEAMS_CHILL_PUBLISHED,
 *                             labels, n, NULL);
 *   }
 *   if (rc != SEAMS_OK) {
 *     fprintf(stderr, "%s: %s\n", seams_status_string(rc), seams_last_error_message());
 *   }
 *   seams_system_destroy(sys);
 *   return (int)rc;
 * }
 * @endcode
 *
 * The version macros live in seams_version.h, which this header
 * includes: SEAMS_VERSION_MAJOR, SEAMS_VERSION_MINOR,
 * SEAMS_VERSION_PATCH, SEAMS_VERSION_STRING, SEAMS_VERSION_NUMBER.
 *
 * SPDX-License-Identifier: MIT
 * Copyright (c) d-SEAMS developers
 */

/** @brief ABI generation of this header. Bumped only on a breaking change. */
#define SEAMS_ABI 1
#include "seams_version.h"

#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @def SEAMS_API
 * @brief Marks a symbol exported by libseams.
 *
 * The library is compiled with hidden visibility and SEAMS_BUILDING
 * defined, so only declarations carrying this macro are exported.
 * Windows consumers get dllimport. Define SEAMS_STATIC when linking a
 * static archive.
 */
#if defined(SEAMS_STATIC)
#define SEAMS_API
#elif defined(_WIN32) || defined(__CYGWIN__)
#if defined(SEAMS_BUILDING)
#define SEAMS_API __declspec(dllexport)
#else
#define SEAMS_API __declspec(dllimport)
#endif
#elif defined(__GNUC__) || defined(__clang__)
#define SEAMS_API __attribute__((visibility("default")))
#else
#define SEAMS_API
#endif

/**
 * @def SEAMS_DEPRECATED
 * @brief Marks an entry point kept for ABI 1 callers. @p msg names the replacement.
 */
#if defined(__GNUC__) || defined(__clang__)
#define SEAMS_DEPRECATED(msg) __attribute__((deprecated(msg)))
#elif defined(_MSC_VER)
#define SEAMS_DEPRECATED(msg) __declspec(deprecated(msg))
#else
#define SEAMS_DEPRECATED(msg)
#endif

/**
 * @defgroup capi_types Types
 * @brief Status codes, label codes, and opaque handles.
 * @{
 */

/** @brief Status codes. Zero is success. */
typedef enum seams_status {
  SEAMS_OK = 0,            /**< The call succeeded. */
  SEAMS_ERR_NULL = 1,      /**< A required pointer is null. */
  SEAMS_ERR_RANGE = 2,     /**< A number or a buffer length is outside the accepted range. */
  SEAMS_ERR_BOX = 3,       /**< The cell is missing, non-finite, or has no volume. */
  SEAMS_ERR_ARG = 4,       /**< A caller argument is the wrong shape, or a coordinate is not finite. */
  SEAMS_ERR_STATE = 5,     /**< The handle has no frame yet. */
  SEAMS_ERR_NOMEM = 6,     /**< An allocation failed. */
  SEAMS_ERR_INTERNAL = 7   /**< The library hit an internal failure. */
} seams_status;

/** @brief CHILL and CHILL+ atom labels. Matches molSys::atom_state_type. */
typedef enum seams_label {
  SEAMS_LABEL_CUBIC = 0,          /**< Cubic ice (Ic). */
  SEAMS_LABEL_HEXAGONAL = 1,      /**< Hexagonal ice (Ih). */
  SEAMS_LABEL_WATER = 2,          /**< Liquid water. */
  SEAMS_LABEL_INTERFACIAL = 3,    /**< Interfacial ice. */
  SEAMS_LABEL_CLATHRATE = 4,      /**< Clathrate. */
  SEAMS_LABEL_INTERCLATHRATE = 5, /**< Interfacial clathrate. */
  SEAMS_LABEL_UNCLASSIFIED = 6,   /**< Not classified, or dropped by the species filter. */
  SEAMS_LABEL_RECUBIC = 7,        /**< Cubic after reclassification. */
  SEAMS_LABEL_REHEX = 8           /**< Hexagonal after reclassification. */
} seams_label;

/** @brief Ring topology. Size classes are used when a six-ring is not a cage face. */
typedef enum seams_ring_symmetry {
  SEAMS_RING_UNCLASSIFIED = 0,   /**< No topology assigned. */
  SEAMS_RING_DDC = 1,            /**< Face of a double-diamond cage. */
  SEAMS_RING_HC_BASAL = 2,       /**< Basal face of a hexagonal cage. */
  SEAMS_RING_HC_PRISMATIC = 3,   /**< Prismatic face of a hexagonal cage. */
  SEAMS_RING_BOTH_BASAL = 4,     /**< Shared by both cage kinds, basal. */
  SEAMS_RING_BOTH_PRISMATIC = 5, /**< Shared by both cage kinds, prismatic. */
  SEAMS_RING_PRISM = 6,          /**< Face of a perfect prism. */
  SEAMS_RING_DEFORMED_PRISM = 7, /**< Face of a deformed prism. */
  SEAMS_RING_MIXED_PRISM = 8,    /**< Face of a mixed prism. */
  SEAMS_RING_PENTAGON = 9,       /**< Five-ring outside any cage. */
  SEAMS_RING_HEXAGON = 10,       /**< Six-ring outside any cage. */
  SEAMS_RING_OTHER = 11          /**< Any other size outside any cage. */
} seams_ring_symmetry;

/**
 * @brief CHILL+ definition.
 *
 * PUBLISHED averages q_lm over the four closest neighbours (Nguyen and
 * Molinero, equation 1) and uses inclusive edges.
 * OVITO averages q_lm over every neighbour inside the cutoff, the Chill+
 * modifier's 3.5 A shell, with that modifier's open edges.
 */
typedef enum seams_chill_rule {
  SEAMS_CHILL_PUBLISHED = 0, /**< Nguyen and Molinero, four closest neighbours. */
  SEAMS_CHILL_OVITO = 1      /**< OVITO Chill+ modifier, every neighbour in the cutoff. */
} seams_chill_rule;

/** @brief Neighbour graph. Seeded is a cage graph: mutual seeds, union completions. */
typedef enum seams_graph {
  SEAMS_GRAPH_CUTOFF = 0, /**< Every neighbour inside the cutoff. */
  SEAMS_GRAPH_MUTUAL = 1, /**< Mutual four nearest inside cutoff + 1.5. */
  SEAMS_GRAPH_SEEDED = 2  /**< Mutual seeds with union completions, for cages only. */
} seams_graph;

/**
 * @brief Per-atom cage membership from the hexagonal and double-diamond search.
 *
 * The four codes are disjoint.
 */
typedef enum seams_cage_label {
  SEAMS_CAGE_WATER = 0, /**< In neither cage. */
  SEAMS_CAGE_HEX = 1,   /**< Hexagonal cage only. */
  SEAMS_CAGE_CUBIC = 2, /**< Double-diamond cage only. */
  SEAMS_CAGE_BOTH = 3   /**< Both cages. */
} seams_cage_label;

/** @brief Site kind. Values match the engine's site table. */
typedef enum seams_site_kind {
  SEAMS_SITE_UNSPECIFIED = 0, /**< Not named in the site table. */
  SEAMS_SITE_CATION = 1,      /**< Cation head. */
  SEAMS_SITE_ANION = 2,       /**< Anion. */
  SEAMS_SITE_TAIL = 3,        /**< Alkyl tail. */
  SEAMS_SITE_DONOR_H = 4,     /**< Hydrogen-bond donor hydrogen. */
  SEAMS_SITE_ACCEPTOR = 5,    /**< Hydrogen-bond acceptor. */
  SEAMS_SITE_POLAR = 6,       /**< Other polar site. */
  SEAMS_SITE_APOLAR = 7,      /**< Other apolar site. */
  SEAMS_SITE_WATER_O = 8,     /**< Water oxygen. */
  SEAMS_SITE_WATER_H = 9,     /**< Water hydrogen. */
  SEAMS_SITE_SOLVENT = 10     /**< Other solvent. */
} seams_site_kind;

/** @brief One frame. Create with seams_system_create(), free with seams_system_destroy(). */
typedef struct seams_system seams_system;
/** @brief One neighbour list. Free with seams_neighbours_destroy(). */
typedef struct seams_neighbours seams_neighbours;
/** @brief One primitive-ring set. Free with seams_rings_destroy(). */
typedef struct seams_rings seams_rings;
/** @brief Hexagonal and double-diamond membership. Free with seams_cages_destroy(). */
typedef struct seams_cage_set seams_cage_set;
/** @brief One planar ice/brine front measurement. Free with seams_brine_destroy(). */
typedef struct seams_brine seams_brine;
/** @brief Topological fingerprint of one bonded frame. Free with seams_fingerprint_destroy(). */
typedef struct seams_topo seams_topo;

/** @} */

/**
 * @defgroup capi_meta Library metadata and errors
 * @brief Version, names, error strings, and the thread team.
 * @{
 */

/**
 * @brief ABI number compiled into this library.
 * @return SEAMS_ABI, which is 1.
 * @par Thread safety
 * Pure.
 * @par Example
 * @code
 * if (seams_abi() != SEAMS_ABI) { return 1; }
 * @endcode
 */
SEAMS_API int seams_abi(void);

/**
 * @brief Packed release number.
 * @return SEAMS_VERSION_NUMBER, (major << 16) | (minor << 8) | patch.
 * @par Thread safety
 * Pure.
 * @par Example
 * @code
 * if (seams_version() < ((2 << 16) | (14 << 8))) { return 1; }
 * @endcode
 */
SEAMS_API int seams_version(void);

/**
 * @brief Major version of the library at run time.
 * @return SEAMS_VERSION_MAJOR of the build.
 * @par Thread safety
 * Pure.
 * @par Example
 * @code
 * int major = seams_version_major();
 * @endcode
 */
SEAMS_API int seams_version_major(void);

/**
 * @brief Minor version of the library at run time.
 * @return SEAMS_VERSION_MINOR of the build.
 * @par Thread safety
 * Pure.
 * @par Example
 * @code
 * int minor = seams_version_minor();
 * @endcode
 */
SEAMS_API int seams_version_minor(void);

/**
 * @brief Patch version of the library at run time.
 * @return SEAMS_VERSION_PATCH of the build.
 * @par Thread safety
 * Pure.
 * @par Example
 * @code
 * int patch = seams_version_patch();
 * @endcode
 */
SEAMS_API int seams_version_patch(void);

/**
 * @brief Version text, for example "2.14.0".
 * @return Pointer to a static string. Never NULL.
 * @par Ownership
 * The caller does not free it.
 * @par Thread safety
 * Pure.
 * @par Example
 * @code
 * printf("d-SEAMS %s\n", seams_version_string());
 * @endcode
 */
SEAMS_API const char *seams_version_string(void);

/**
 * @brief Short name of a status code.
 * @param status A ::seams_status value.
 * @return Pointer to a static string. Unknown codes return "unknown status".
 * @par Ownership
 * The caller does not free it.
 * @par Thread safety
 * Pure.
 * @par Example
 * @code
 * fprintf(stderr, "%s\n", seams_status_string(SEAMS_ERR_BOX));
 * @endcode
 */
SEAMS_API const char *seams_status_string(seams_status status);

/**
 * @brief Thread-local detail of the last failing call on this thread.
 *
 * ABI 1 spelling of seams_last_error_message(). Both return the same
 * string.
 *
 * @return Pointer to a thread-local string. Never NULL.
 * @par Ownership
 * Owned by the library. Valid until the next call on this thread that
 * sets or clears the error.
 * @par Thread safety
 * The storage is thread-local.
 * @par Example
 * @code
 * if (seams_system_set_atoms(sys, xyz, n, box9, 9, NULL, NULL) != SEAMS_OK) {
 *   fprintf(stderr, "%s\n", seams_last_error());
 * }
 * @endcode
 */
SEAMS_API const char *seams_last_error(void);

/**
 * @brief Detail string stored on a handle.
 * @param handle A ::seams_system, ::seams_neighbours, ::seams_rings,
 *        ::seams_cage_set, ::seams_topo, or ::seams_brine, or NULL.
 * @return The handle's string, or the thread-local string when handle
 *         is NULL. Never NULL.
 *
 * A call taking the handle non-const clears the string on entry. A call
 * taking it const writes the string only when it fails, so after a
 * successful const call the string still holds the last failure.
 * @par Ownership
 * Owned by the handle. Valid until the next failing call on that
 * handle, the next non-const call on it, or destroy.
 * @par Thread safety
 * Reading a handle requires that no other thread is inside a call on
 * that same handle.
 * @par Example
 * @code
 * if (seams_chill_plus_n(sys, 3.5, 0, SEAMS_GRAPH_CUTOFF, SEAMS_CHILL_PUBLISHED,
 *                        labels, n, NULL) != SEAMS_OK) {
 *   fprintf(stderr, "%s\n", seams_handle_error(sys));
 * }
 * @endcode
 */
SEAMS_API const char *seams_handle_error(const void *handle);

/**
 * @brief Name of a CHILL label.
 * @param label A ::seams_label value.
 * @return Static string. Unknown codes return "unknown".
 * @par Ownership
 * The caller does not free it.
 * @par Thread safety
 * Pure.
 * @par Example
 * @code
 * printf("%s\n", seams_label_name(SEAMS_LABEL_HEXAGONAL));
 * @endcode
 */
SEAMS_API const char *seams_label_name(int label);

/**
 * @brief Name of a ring-symmetry code.
 * @param symmetry A ::seams_ring_symmetry value.
 * @return Static string. Unknown codes return "unknown".
 * @par Ownership
 * The caller does not free it.
 * @par Thread safety
 * Pure.
 * @par Example
 * @code
 * printf("%s\n", seams_ring_symmetry_name(SEAMS_RING_DDC));
 * @endcode
 */
SEAMS_API const char *seams_ring_symmetry_name(int symmetry);

/**
 * @brief Set the OpenMP team used by later analysis calls.
 * @param n Positive width, or 0 to restore the width captured at first
 *        use. A negative width is rejected.
 * @return ::SEAMS_OK, or ::SEAMS_ERR_RANGE when n is negative,
 *         ::SEAMS_ERR_NOMEM, ::SEAMS_ERR_INTERNAL.
 * @par Thread safety
 * Process-global. Do not call it concurrently with analysis. The stored
 * width is then read by every handle.
 * @par Example
 * @code
 * seams_set_num_threads(1);
 * @endcode
 */
SEAMS_API seams_status seams_set_num_threads(int n);

/**
 * @brief Team width the next analysis call will request.
 * @return A positive count. After seams_set_num_threads(0) this is the
 *         saved default, which is omp_get_max_threads() when OpenMP is
 *         linked and 1 otherwise.
 * @par Thread safety
 * Observes the process-global width. Do not call it concurrently with
 * seams_set_num_threads().
 * @par Example
 * @code
 * int width = seams_get_num_threads();
 * @endcode
 */
SEAMS_API int seams_get_num_threads(void);

/** @} */

/**
 * @defgroup capi_frame Frames
 * @brief Create, fill, read, and free a ::seams_system.
 * @{
 */

/**
 * @brief Allocate an empty frame.
 * @return A new handle, or NULL on allocation failure.
 *         seams_last_error_message() then explains the failure.
 * @par Ownership
 * The caller owns the handle and passes it to seams_system_destroy().
 * @par Thread safety
 * Safe. The new handle is not shared yet.
 * @par Example
 * @code
 * seams_system *sys = seams_system_create();
 * if (sys == NULL) { return SEAMS_ERR_NOMEM; }
 * seams_system_destroy(sys);
 * @endcode
 */
SEAMS_API seams_system *seams_system_create(void);

/**
 * @brief Free a frame.
 * @param sys Handle from seams_system_create(), or NULL. NULL is ignored.
 * @par Ownership
 * The handle and the frame it stored are released. Neighbour, ring,
 * cage, and brine handles built from it stay valid: they hold their
 * own copies.
 * @par Thread safety
 * The caller must not use sys after this returns, and must not call it
 * concurrently with any use of sys.
 * @par Example
 * @code
 * seams_system_destroy(sys);
 * sys = NULL;
 * @endcode
 */
SEAMS_API void seams_system_destroy(seams_system *sys);

/**
 * @brief Store points and a general cell. Lengths are size_t.
 *
 * On failure the previously stored frame, if any, is kept.
 *
 * A reload on the same handle reuses its per-atom storage. When @p ids
 * is NULL or 1..n and @p n matches the previous load, the id table is
 * kept as well, so reloading a trajectory frame by frame costs one pass
 * over the input after the finiteness check. seams_system_set_positions()
 * also skips the type and id pass.
 *
 * @param sys Frame handle.
 * @param xyz n times 3 Cartesian coordinates, Angstrom. Every value must
 *        be finite.
 * @param n Atom count. Zero is rejected. A count that does not fit in
 *        the frame's integer width is rejected.
 * @param box9 Cell numbers: ox, oy, oz, lx, ly, lz, xy, xz, yz.
 * @param box_len Must be 9.
 * @param types Atom types, or NULL to store type 1 for every atom.
 * @param ids Atom ids, or NULL to store 1..n. Ids must be unique.
 * @return ::SEAMS_OK; ::SEAMS_ERR_NULL for a NULL sys, xyz, or box9;
 *         ::SEAMS_ERR_RANGE for n == 0 or box_len != 9; ::SEAMS_ERR_BOX
 *         for a non-finite cell or a non-positive edge; ::SEAMS_ERR_ARG
 *         for a non-finite coordinate or repeated ids; ::SEAMS_ERR_NOMEM;
 *         ::SEAMS_ERR_INTERNAL.
 * @par Handle access
 * Non-const: replaces the stored frame. On failure the old frame stays.
 * @par Ownership
 * The library copies the arrays. The caller keeps them.
 * @par Thread safety
 * Exclusive use of sys.
 * @par Example
 * @code
 * const double xyz[6] = {0, 0, 0, 1.5, 0, 0};
 * const double box9[9] = {0, 0, 0, 10, 10, 10, 0, 0, 0};
 * seams_status rc = seams_system_set_atoms(sys, xyz, 2, box9, 9, NULL, NULL);
 * @endcode
 */
SEAMS_API seams_status seams_system_set_atoms(seams_system *sys, const double *xyz,
                                              size_t n, const double *box9,
                                              size_t box_len, const int *types,
                                              const int *ids);

/**
 * @brief Store points and a general LAMMPS cell, with an int count.
 *
 * ABI 1 entry. New callers use seams_system_set_atoms().
 *
 * @param sys Frame handle.
 * @param xyz n times 3 coordinates, Angstrom.
 * @param n Atom count. A non-positive count is rejected.
 * @param box9 Nine cell numbers: ox, oy, oz, lx, ly, lz, xy, xz, yz.
 * @param types Types, or NULL for type 1.
 * @param ids Ids, or NULL for 1..n.
 * @return Same codes as seams_system_set_atoms().
 * @par Handle access
 * Non-const: replaces the stored frame. On failure the old frame stays.
 * @par Ownership
 * The library copies the arrays.
 * @par Thread safety
 * Exclusive use of sys.
 * @par Example
 * @code
 * seams_status rc = seams_system_set_frame(sys, xyz, 2, box9, NULL, NULL);
 * @endcode
 */
SEAMS_API seams_status seams_system_set_frame(seams_system *sys, const double *xyz,
                                              int n, const double *box9,
                                              const int *types, const int *ids);

/**
 * @brief Store points and an orthorhombic cell.
 * @param sys Frame handle.
 * @param xyz n times 3 coordinates, Angstrom.
 * @param n Atom count.
 * @param ox Origin x, Angstrom.
 * @param oy Origin y, Angstrom.
 * @param oz Origin z, Angstrom.
 * @param lx Edge length x, Angstrom.
 * @param ly Edge length y, Angstrom.
 * @param lz Edge length z, Angstrom.
 * @param types Types, or NULL.
 * @param ids Ids, or NULL.
 * @return Same codes as seams_system_set_frame(). The tilt is zero.
 * @par Handle access
 * Non-const: replaces the stored frame. On failure the old frame stays.
 * @par Ownership
 * The library copies the arrays.
 * @par Thread safety
 * Exclusive use of sys.
 * @par Example
 * @code
 * seams_status rc = seams_system_set_ortho(sys, xyz, 2, 0, 0, 0, 10, 10, 10, NULL, NULL);
 * @endcode
 */
SEAMS_API seams_status seams_system_set_ortho(seams_system *sys, const double *xyz,
                                              int n, double ox, double oy, double oz,
                                              double lx, double ly, double lz,
                                              const int *types, const int *ids);

/**
 * @brief Read one frame from a trajectory file.
 * @param sys Frame handle. On failure the previously stored frame is kept.
 * @param path Filesystem path. CON, LAMMPS dump, and the formats the
 *        linked readers accept.
 * @param frame Frame index, 1-based, matching the file order.
 * @param type_filter Negative keeps every atom. Zero guesses the oxygen
 *        species. A positive value keeps that species number.
 * @return ::SEAMS_OK, or ::SEAMS_ERR_NULL, ::SEAMS_ERR_RANGE (frame index
 *         or empty frame), ::SEAMS_ERR_ARG (unreadable file),
 *         ::SEAMS_ERR_NOMEM, ::SEAMS_ERR_INTERNAL.
 * @par Handle access
 * Non-const: replaces the stored frame. On failure the old frame stays.
 * @par Ownership
 * The library copies the frame out of the file.
 * @par Thread safety
 * Exclusive use of sys. Concurrent reads of distinct paths on distinct
 * handles are safe.
 * @par Example
 * @code
 * seams_status rc = seams_system_read(sys, "ice.lammpstrj", 1, -1);
 * @endcode
 */
SEAMS_API seams_status seams_system_read(seams_system *sys, const char *path,
                                         int frame, int type_filter);

/**
 * @brief Number of frames in a trajectory path.
 * @param path Filesystem path.
 * @param nframes Receives the count. Must be non-NULL.
 * @return ::SEAMS_OK, or ::SEAMS_ERR_NULL, ::SEAMS_ERR_RANGE,
 *         ::SEAMS_ERR_ARG, ::SEAMS_ERR_INTERNAL. *nframes is unchanged on
 *         failure.
 * @par Ownership
 * path is read during the call only. nframes is caller storage that the call
 * fills. The library keeps no pointer to any argument.
 * @par Thread safety
 * Safe for distinct paths.
 * @par Example
 * @code
 * size_t nframes = 0;
 * seams_status rc = seams_path_nframes("ice.lammpstrj", &nframes);
 * @endcode
 */
SEAMS_API seams_status seams_path_nframes(const char *path, size_t *nframes);

/**
 * @brief Number of points in the frame, as an int.
 *
 * ABI 1 query. It cannot report an error; seams_system_natoms() can.
 *
 * @param sys Frame handle, or NULL.
 * @return The atom count, or 0 when sys is NULL or holds no frame.
 * @par Ownership
 * The caller owns sys; the call borrows it. The library keeps no pointer to
 * any argument.
 * @par Thread safety
 * Read-only: may run alongside other const calls on sys, and not
 * alongside a non-const call on sys.
 * @par Example
 * @code
 * int n = seams_system_count(sys);
 * @endcode
 */
SEAMS_API int seams_system_count(const seams_system *sys);

/** @} */

/**
 * @defgroup capi_neigh Neighbour lists
 * @brief Cutoff, mutual, and k-nearest graphs as CSR.
 * @{
 */

/**
 * @brief Build a cutoff or mutual neighbour list.
 * @param sys Frame handle that already holds a frame.
 * @param cutoff Positive finite cutoff, Angstrom. The mutual graph keeps
 *        the four nearest neighbours inside cutoff + 1.5.
 * @param type Species to keep. Zero and negative keep every atom.
 * @param graph ::SEAMS_GRAPH_CUTOFF or ::SEAMS_GRAPH_MUTUAL.
 *        ::SEAMS_GRAPH_SEEDED is rejected here.
 * @param out Receives a new handle. On failure *out is NULL.
 * @return ::SEAMS_OK, or ::SEAMS_ERR_NULL, ::SEAMS_ERR_STATE,
 *         ::SEAMS_ERR_RANGE, ::SEAMS_ERR_ARG, ::SEAMS_ERR_NOMEM,
 *         ::SEAMS_ERR_INTERNAL.
 * @par Handle access
 * Non-const only because the call clears the handle's error slot on
 * entry, so seams_handle_error() reports this call. It reads the frame
 * and keeps no scratch on the handle.
 * @par Ownership
 * The caller frees *out with seams_neighbours_destroy(). The list is a
 * copy and does not borrow sys.
 * @par Thread safety
 * Exclusive use of sys. The new handle is private to the caller.
 * @par Example
 * @code
 * seams_neighbours *nb = NULL;
 * if (seams_neighbours_build(sys, 3.5, 0, SEAMS_GRAPH_CUTOFF, &nb) == SEAMS_OK) {
 *   seams_neighbours_destroy(nb);
 * }
 * @endcode
 */
SEAMS_API seams_status seams_neighbours_build(seams_system *sys, double cutoff,
                                              int type, int graph,
                                              seams_neighbours **out);

/**
 * @brief Free a neighbour list.
 * @param nb Handle from seams_neighbours_build() or seams_neighbours_knn(),
 *        or NULL. NULL is ignored.
 * @par Ownership
 * Releases the handle and its CSR arrays.
 * @par Thread safety
 * The caller must not use nb after this returns.
 * @par Example
 * @code
 * seams_neighbours_destroy(nb);
 * @endcode
 */
SEAMS_API void seams_neighbours_destroy(seams_neighbours *nb);

/**
 * @brief Number of neighbour rows, as an int.
 *
 * ABI 1 query. seams_neighbours_nrows() reports errors.
 *
 * @param nb Neighbour handle, or NULL.
 * @return The row count, or 0 when nb is NULL.
 * @par Ownership
 * The caller owns nb; the call borrows it. The library keeps no pointer to
 * any argument.
 * @par Thread safety
 * Read-only: may run alongside other const calls on nb, and not
 * alongside a non-const call on nb.
 * @par Example
 * @code
 * int rows = seams_neighbours_count(nb);
 * @endcode
 */
SEAMS_API int seams_neighbours_count(const seams_neighbours *nb);

/**
 * @brief Borrow the neighbour list as CSR.
 *
 * Offsets has n + 1 entries. Columns omit the self entry.
 *
 * @param nb Neighbour handle.
 * @param offsets Receives the offset pointer.
 * @param n_offsets Receives the offset count, n + 1.
 * @param cols Receives the column pointer. May be NULL when the list
 *        has no pairs; *nnz is then 0.
 * @param nnz Receives the column count.
 * @return ::SEAMS_OK, or ::SEAMS_ERR_NULL.
 * @par Ownership
 * The pointers are owned by nb and remain valid until
 * seams_neighbours_destroy(nb). Use seams_neighbours_copy() for caller
 * buffers.
 * @par Thread safety
 * Exclusive use of nb. The returned pointers must not be used
 * concurrently with a call that frees nb.
 * @par Example
 * @code
 * const int *off = NULL, *cols = NULL;
 * int n_off = 0, nnz = 0;
 * seams_neighbours_csr(nb, &off, &n_off, &cols, &nnz);
 * @endcode
 */
SEAMS_DEPRECATED("use seams_neighbours_copy; these pointers are owned by the handle")
SEAMS_API seams_status seams_neighbours_csr(const seams_neighbours *nb,
                                            const int **offsets, int *n_offsets,
                                            const int **cols, int *nnz);

/**
 * @brief Copy the neighbour CSR into caller buffers.
 * @param nb Neighbour handle.
 * @param offsets Caller buffer of offsets_cap ints, or NULL to query.
 * @param offsets_cap Capacity of offsets. Ignored when offsets is NULL.
 * @param n_offsets Receives the required offset count (n + 1).
 * @param cols Caller buffer of cols_cap ints, or NULL to query.
 * @param cols_cap Capacity of cols. Ignored when cols is NULL.
 * @param nnz Receives the required column count.
 * @return ::SEAMS_OK when every non-NULL buffer was filled, or when both
 *         buffers were NULL and the lengths were written.
 *         ::SEAMS_ERR_RANGE when a non-NULL buffer is short.
 *         ::SEAMS_ERR_NULL when nb, n_offsets, or nnz is NULL.
 * @par Ownership
 * The caller owns the buffers. The handle is not modified.
 * @par Thread safety
 * Read-only: may run alongside other const calls on nb, and not
 * alongside a non-const call on nb.
 * @par Example
 * @code
 * size_t n_off = 0, nnz = 0;
 * seams_neighbours_copy(nb, NULL, 0, &n_off, NULL, 0, &nnz);
 * int *off = malloc(n_off * sizeof *off), *cols = malloc(nnz * sizeof *cols);
 * seams_neighbours_copy(nb, off, n_off, &n_off, cols, nnz, &nnz);
 * @endcode
 */
SEAMS_API seams_status seams_neighbours_copy(const seams_neighbours *nb,
                                             int *offsets, size_t offsets_cap,
                                             size_t *n_offsets, int *cols,
                                             size_t cols_cap, size_t *nnz);

/** @} */

/**
 * @defgroup capi_order Order parameters
 * @brief CHILL, CHILL+, Steinhardt q_l, and the ice shell label.
 * @{
 */

/**
 * @brief Original CHILL labels, published or OVITO edges.
 *
 * Published CHILL is Moore's four-neighbour rule with inclusive c_ij
 * edges. OVITO mode uses the same thresholds and the open edges.
 * Labels match ::seams_label.
 *
 * @param sys Frame handle.
 * @param cutoff Positive finite cutoff, Angstrom.
 * @param type Species filter. Zero and negative keep every atom.
 * @param graph ::SEAMS_GRAPH_CUTOFF or ::SEAMS_GRAPH_MUTUAL.
 * @param rule ::SEAMS_CHILL_PUBLISHED or ::SEAMS_CHILL_OVITO.
 * @param labels Caller buffer of labels_len ints.
 * @param labels_len Must be at least the frame count.
 * @param counts9 Nine ints that receive the population of each label,
 *        or NULL to skip the histogram.
 * @return ::SEAMS_OK, or ::SEAMS_ERR_NULL, ::SEAMS_ERR_STATE,
 *         ::SEAMS_ERR_RANGE, ::SEAMS_ERR_ARG, ::SEAMS_ERR_NOMEM,
 *         ::SEAMS_ERR_INTERNAL.
 * @par Handle access
 * Non-const: writes each atom's CHILL+ label into the stored frame and
 * clears the handle's error slot on entry. Coordinates, cell, types, and
 * ids are not changed.
 * @par Ownership
 * The caller owns labels and counts9.
 * @par Thread safety
 * Exclusive use of sys.
 * @par Example
 * @code
 * int counts[9];
 * seams_chill(sys, 3.5, 0, SEAMS_GRAPH_CUTOFF, SEAMS_CHILL_PUBLISHED, labels, n, counts);
 * @endcode
 */
SEAMS_API seams_status seams_chill(seams_system *sys, double cutoff, int type,
                                   int graph, int rule, int *labels,
                                   size_t labels_len, int counts9[9]);

/**
 * @brief CHILL+ labels, without a buffer length.
 *
 * ABI 1 entry. labels must hold one int per atom in the frame. New
 * callers use seams_chill_plus_n(), which checks the length.
 *
 * @param sys Frame handle.
 * @param cutoff Positive finite cutoff, Angstrom.
 * @param type Species filter. Zero and negative keep every atom.
 * @param graph ::SEAMS_GRAPH_CUTOFF or ::SEAMS_GRAPH_MUTUAL.
 * @param rule ::SEAMS_CHILL_PUBLISHED or ::SEAMS_CHILL_OVITO.
 * @param labels Caller buffer of one label per atom.
 * @param counts9 Nine ints, or NULL.
 * @return Same codes as seams_chill_plus_n().
 * @par Handle access
 * Non-const: writes each atom's CHILL+ label into the stored frame and
 * clears the handle's error slot on entry. Coordinates, cell, types, and
 * ids are not changed.
 * @par Ownership
 * The caller owns labels and counts9.
 * @par Thread safety
 * Exclusive use of sys.
 * @par Example
 * @code
 * seams_chill_plus(sys, 3.5, 0, SEAMS_GRAPH_CUTOFF, SEAMS_CHILL_PUBLISHED, labels, NULL);
 * @endcode
 */
SEAMS_API seams_status seams_chill_plus(seams_system *sys, double cutoff, int type,
                                        int graph, int rule, int *labels,
                                        int counts9[9]);

/**
 * @brief Local Steinhardt q_l, and the Lechner-Dellago average on request.
 *
 * order_l is 3, 4, 6, 8, or 12. Either array may be NULL. A non-NULL
 * array has one double per atom and its length must be at least the
 * frame count. The average is the neighbour mean of q_lm. The cutoff
 * graph reads the cutoff shell. The mutual graph reads the four nearest
 * neighbours inside cutoff + 1.5.
 *
 * @param sys Frame handle.
 * @param cutoff Positive finite cutoff, Angstrom.
 * @param type Species filter. Zero and negative keep every atom.
 * @param order_l Spherical-harmonic degree.
 * @param graph ::SEAMS_GRAPH_CUTOFF or ::SEAMS_GRAPH_MUTUAL.
 * @param ql Local q_l, or NULL.
 * @param ql_len Length of ql. Ignored when ql is NULL.
 * @param qlbar Averaged q_l, or NULL to skip the average.
 * @param qlbar_len Length of qlbar. Ignored when qlbar is NULL.
 * @return ::SEAMS_OK, or ::SEAMS_ERR_NULL, ::SEAMS_ERR_STATE,
 *         ::SEAMS_ERR_RANGE, ::SEAMS_ERR_ARG, ::SEAMS_ERR_NOMEM,
 *         ::SEAMS_ERR_INTERNAL.
 * @par Handle access
 * Non-const only because the call clears the handle's error slot on
 * entry; it reads the frame and nothing else on the handle. The q_lm work
 * arrays belong to the calling thread, so calls on different handles may
 * run at the same time.
 * @par Ownership
 * The caller owns both arrays.
 * @par Thread safety
 * Exclusive use of sys. Any q call on another handle may run beside it.
 * @par Example
 * @code
 * seams_ql(sys, 3.5, 0, 6, SEAMS_GRAPH_CUTOFF, q6, n, q6bar, n);
 * @endcode
 */
SEAMS_API seams_status seams_ql(seams_system *sys, double cutoff, int type,
                                int order_l, int graph, double *ql, size_t ql_len,
                                double *qlbar, size_t qlbar_len);

/**
 * @brief Local q4 and q6 on the cutoff shell, without buffer lengths.
 *
 * ABI 1 pair. Each non-NULL array has one double per atom. New callers
 * use seams_q4q6_n(), which checks the lengths.
 *
 * @param sys Frame handle.
 * @param cutoff Positive finite cutoff, Angstrom.
 * @param type Species filter.
 * @param q4 Local q4, or NULL.
 * @param q6 Local q6, or NULL. q4 and q6 may not both be NULL.
 * @return ::SEAMS_OK, or ::SEAMS_ERR_NULL, ::SEAMS_ERR_STATE,
 *         ::SEAMS_ERR_RANGE, ::SEAMS_ERR_NOMEM, ::SEAMS_ERR_INTERNAL.
 * @par Handle access
 * Non-const only because the call clears the handle's error slot on
 * entry; it reads the frame and nothing else on the handle. The q_lm work
 * arrays belong to the calling thread, so calls on different handles may
 * run at the same time.
 * @par Ownership
 * The caller owns q4 and q6.
 * @par Thread safety
 * Exclusive use of sys. Any q call on another handle may run beside it.
 * @par Example
 * @code
 * seams_q4q6(sys, 3.5, 0, q4, q6);
 * @endcode
 */
SEAMS_API seams_status seams_q4q6(seams_system *sys, double cutoff, int type,
                                  double *q4, double *q6);

/**
 * @brief Lechner-Dellago q6-bar and the ice count above a cut.
 *
 * The average is the same neighbour mean of q_6m as seams_ql() of
 * order 6 on the cutoff graph. n_ice counts finite values strictly
 * above ice_cut. The nucleation census uses 0.36. mean and median
 * ignore non-finite values. An empty finite set stores 0.
 *
 * @param sys Frame handle.
 * @param cutoff Positive finite cutoff, Angstrom.
 * @param type Species filter.
 * @param ice_cut Finite threshold on q6-bar.
 * @param q6bar Per-atom average, or NULL to keep only the summary.
 * @param q6bar_len Length of q6bar. Ignored when q6bar is NULL.
 * @param mean Receives the mean, or NULL.
 * @param median Receives the median, or NULL.
 * @param n_ice Receives the count above ice_cut, or NULL.
 * @return ::SEAMS_OK, or ::SEAMS_ERR_NULL, ::SEAMS_ERR_STATE,
 *         ::SEAMS_ERR_RANGE, ::SEAMS_ERR_NOMEM, ::SEAMS_ERR_INTERNAL. At
 *         least one of q6bar, mean, median, and n_ice must be non-NULL.
 * @par Handle access
 * Non-const only because the call clears the handle's error slot on
 * entry; it reads the frame and nothing else on the handle. The q_lm work
 * arrays belong to the calling thread, so calls on different handles may
 * run at the same time.
 * @par Ownership
 * The caller owns q6bar.
 * @par Thread safety
 * Exclusive use of sys. Any q call on another handle may run beside it.
 * @par Example
 * @code
 * double mean = 0, median = 0;
 * int n_ice = 0;
 * seams_q6bar(sys, 3.5, 0, 0.36, NULL, 0, &mean, &median, &n_ice);
 * @endcode
 */
SEAMS_API seams_status seams_q6bar(seams_system *sys, double cutoff, int type,
                                   double ice_cut, double *q6bar, size_t q6bar_len,
                                   double *mean, double *median, int *n_ice);

/** @} */

/**
 * @defgroup capi_ring Rings
 * @brief Primitive rings (Franzblau) and their labels.
 * @{
 */

/**
 * @brief Primitive rings (Franzblau) up to a maximum size.
 * @param sys Frame handle.
 * @param cutoff Positive finite cutoff, Angstrom.
 * @param type Species filter.
 * @param graph ::SEAMS_GRAPH_CUTOFF or ::SEAMS_GRAPH_MUTUAL.
 * @param max_depth Largest ring size, from 3 through 12.
 * @param out Receives a new handle. On failure *out is NULL.
 * @return ::SEAMS_OK, or ::SEAMS_ERR_NULL, ::SEAMS_ERR_STATE,
 *         ::SEAMS_ERR_RANGE, ::SEAMS_ERR_ARG, ::SEAMS_ERR_NOMEM,
 *         ::SEAMS_ERR_INTERNAL.
 * @par Handle access
 * Non-const only because the call clears the handle's error slot on
 * entry, so seams_handle_error() reports this call. It reads the frame
 * and keeps no scratch on the handle.
 * @par Ownership
 * The caller frees *out with seams_rings_destroy().
 * @par Thread safety
 * Exclusive use of sys.
 * @par Example
 * @code
 * seams_rings *rings = NULL;
 * if (seams_rings_find(sys, 3.5, 0, SEAMS_GRAPH_CUTOFF, 6, &rings) == SEAMS_OK) {
 *   seams_rings_destroy(rings);
 * }
 * @endcode
 */
SEAMS_API seams_status seams_rings_find(seams_system *sys, double cutoff, int type,
                                        int graph, int max_depth, seams_rings **out);

/**
 * @brief Free a ring set.
 * @param rings Handle from seams_rings_find(), or NULL. NULL is ignored.
 * @par Ownership
 * Releases the handle and every ring it holds.
 * @par Thread safety
 * The caller must not use rings after this returns.
 * @par Example
 * @code
 * seams_rings_destroy(rings);
 * @endcode
 */
SEAMS_API void seams_rings_destroy(seams_rings *rings);

/**
 * @brief Number of primitive rings, as an int.
 *
 * ABI 1 query. seams_rings_nrings() reports errors.
 *
 * @param rings Ring handle, or NULL.
 * @return The ring count, or 0 when rings is NULL.
 * @par Ownership
 * The caller owns rings; the call borrows it. The library keeps no pointer to
 * any argument.
 * @par Thread safety
 * Read-only: may run alongside other const calls on rings, and not
 * alongside a non-const call on rings.
 * @par Example
 * @code
 * int nr = seams_rings_count(rings);
 * @endcode
 */
SEAMS_API int seams_rings_count(const seams_rings *rings);

/**
 * @brief Count primitive rings by size, with an int length.
 *
 * ABI 1 entry. New callers use seams_rings_histogram_n().
 *
 * @param rings Ring handle.
 * @param counts Caller buffer. counts[k] receives the number of rings
 *        of size k. Sizes outside the buffer are ignored.
 * @param counts_len Length of counts. Must be positive.
 * @return ::SEAMS_OK, or ::SEAMS_ERR_NULL, ::SEAMS_ERR_RANGE.
 * @par Ownership
 * The caller owns counts.
 * @par Thread safety
 * Read-only: may run alongside other const calls on rings, and not
 * alongside a non-const call on rings.
 * @par Example
 * @code
 * int hist[13];
 * seams_rings_histogram(rings, hist, 13);
 * @endcode
 */
SEAMS_API seams_status seams_rings_histogram(const seams_rings *rings, int *counts,
                                             int counts_len);

/**
 * @brief Size and symmetry of each ring, without buffer lengths.
 *
 * ABI 1 entry. Each non-NULL buffer has length seams_rings_count(rings).
 * New callers use seams_rings_labels_n(), which checks the lengths.
 *
 * @param rings Ring handle.
 * @param sizes One int per ring, or NULL.
 * @param symmetry One ::seams_ring_symmetry code per ring, or NULL.
 *        The two pointers may not both be NULL.
 * @return ::SEAMS_OK, or ::SEAMS_ERR_NULL, ::SEAMS_ERR_NOMEM,
 *         ::SEAMS_ERR_INTERNAL.
 * @par Ownership
 * The caller owns both arrays.
 * @par Thread safety
 * Read-only: may run alongside other const calls on rings, and not
 * alongside a non-const call on rings. The first symmetry request fills
 * a cache on the handle under a lock; later requests only read it.
 * @par Example
 * @code
 * seams_rings_labels(rings, sizes, symmetry);
 * @endcode
 */
SEAMS_API seams_status seams_rings_labels(const seams_rings *rings, int *sizes,
                                          int *symmetry);

/**
 * @brief Atom indices of one ring.
 * @param rings Ring handle.
 * @param index Ring index, from 0 to seams_rings_count(rings) - 1.
 * @param atoms Caller buffer, or NULL to query.
 * @param cap Capacity of atoms. Ignored when atoms is NULL.
 * @param n_out Receives the vertex count.
 * @return ::SEAMS_OK, or ::SEAMS_ERR_NULL, ::SEAMS_ERR_RANGE.
 * @par Ownership
 * The caller owns atoms. Indices are positions in the frame.
 * @par Thread safety
 * Read-only: may run alongside other const calls on rings, and not
 * alongside a non-const call on rings.
 * @par Example
 * @code
 * int atoms[12];
 * size_t len = 0;
 * seams_rings_vertices(rings, 0, atoms, 12, &len);
 * @endcode
 */
SEAMS_API seams_status seams_rings_vertices(const seams_rings *rings, int index,
                                            int *atoms, size_t cap, size_t *n_out);

/**
 * @brief Fraction of atoms that sit on a ring.
 * @param rings Ring handle.
 * @param size Ring size to count. Zero counts atoms that sit on any ring.
 * @param fraction Receives a value in [0, 1].
 * @return ::SEAMS_OK, or ::SEAMS_ERR_NULL, ::SEAMS_ERR_NOMEM.
 * @par Ownership
 * The caller owns rings; the call borrows it. fraction is caller storage that
 * the call fills. The library keeps no pointer to any argument.
 * @par Thread safety
 * Read-only: may run alongside other const calls on rings, and not
 * alongside a non-const call on rings.
 * @par Example
 * @code
 * double f6 = 0;
 * seams_rings_fraction(rings, 6, &f6);
 * @endcode
 */
SEAMS_API seams_status seams_rings_fraction(const seams_rings *rings, int size,
                                            double *fraction);

/** @} */

/**
 * @defgroup capi_cage Cages and ice shells
 * @brief Hexagonal and double-diamond cages, signatures, and shell labels.
 * @{
 */

/**
 * @brief Hexagonal and double-diamond cage membership.
 * @param sys Frame handle.
 * @param cutoff Positive finite cutoff, Angstrom. The seeded graph uses
 *        cutoff + 1.5 as the candidate shell.
 * @param type Species filter.
 * @param graph ::SEAMS_GRAPH_CUTOFF, ::SEAMS_GRAPH_MUTUAL, or
 *        ::SEAMS_GRAPH_SEEDED.
 * @param k Neighbour count for the seeded graph. Zero selects 4.
 * @param complete Non-zero runs ring-adjacent completion on the seeded
 *        labels. Ignored for the cutoff and mutual graphs.
 * @param out Receives a new handle. On failure *out is NULL.
 * @return ::SEAMS_OK, or ::SEAMS_ERR_NULL, ::SEAMS_ERR_STATE,
 *         ::SEAMS_ERR_RANGE, ::SEAMS_ERR_ARG, ::SEAMS_ERR_NOMEM,
 *         ::SEAMS_ERR_INTERNAL.
 * @par Handle access
 * Non-const only because the call clears the handle's error slot on
 * entry, so seams_handle_error() reports this call. It reads the frame
 * and keeps no scratch on the handle.
 * @par Ownership
 * The caller frees *out with seams_cages_destroy().
 * @par Thread safety
 * Exclusive use of sys.
 * @par Example
 * @code
 * seams_cage_set *cages = NULL;
 * if (seams_cages_find(sys, 3.5, 0, SEAMS_GRAPH_CUTOFF, 0, 0, &cages) == SEAMS_OK) {
 *   seams_cages_destroy(cages);
 * }
 * @endcode
 */
SEAMS_API seams_status seams_cages_find(seams_system *sys, double cutoff, int type,
                                        int graph, int k, int complete,
                                        seams_cage_set **out);

/**
 * @brief Free a cage-membership handle.
 * @param cages Handle from seams_cages_find() or seams_cages(), or NULL.
 *        NULL is ignored.
 * @par Ownership
 * Releases the handle and its per-atom labels.
 * @par Thread safety
 * The caller must not use cages after this returns.
 * @par Example
 * @code
 * seams_cages_destroy(cages);
 * @endcode
 */
SEAMS_API void seams_cages_destroy(seams_cage_set *cages);

/**
 * @brief Disjoint cage population counts.
 *
 * The four counts sum to the number of atoms the search kept.
 *
 * @param cages Cage handle.
 * @param n_hex Hexagonal only, or NULL.
 * @param n_cubic Double-diamond only, or NULL.
 * @param n_water Neither cage, or NULL.
 * @param n_both Both cages, or NULL.
 * @return ::SEAMS_OK, or ::SEAMS_ERR_NULL when cages is NULL or every
 *         count pointer is NULL.
 * @par Ownership
 * The caller owns cages; the call borrows it. n_hex, n_cubic, n_water and
 * n_both are caller storage that the call fills. The library keeps no pointer
 * to any argument.
 * @par Thread safety
 * Read-only: may run alongside other const calls on cages, and not
 * alongside a non-const call on cages.
 * @par Example
 * @code
 * int hex = 0, cubic = 0;
 * seams_cages_counts(cages, &hex, &cubic, NULL, NULL);
 * @endcode
 */
SEAMS_API seams_status seams_cages_counts(const seams_cage_set *cages, int *n_hex,
                                          int *n_cubic, int *n_water, int *n_both);

/**
 * @brief Per-atom cage code.
 * @param cages Cage handle.
 * @param per_atom Caller buffer of ::seams_cage_label values, one per atom
 *        in the frame. Atoms the species filter dropped are
 *        ::SEAMS_CAGE_WATER.
 * @param len Length of per_atom. Must be at least the frame count stored
 *        on the handle.
 * @return ::SEAMS_OK, or ::SEAMS_ERR_NULL, ::SEAMS_ERR_RANGE.
 * @par Ownership
 * The caller owns per_atom.
 * @par Thread safety
 * Read-only: may run alongside other const calls on cages, and not
 * alongside a non-const call on cages.
 * @par Example
 * @code
 * seams_cages_labels(cages, per_atom, n);
 * @endcode
 */
SEAMS_API seams_status seams_cages_labels(const seams_cage_set *cages, int *per_atom,
                                          size_t len);

/**
 * @brief Count closed cages whose ring census matches a signature.
 *
 * spec is a census such as "512", "51262", or "51264", or the names
 * "hc" and "ddc". n_cages is the number of closed cages. n_atoms is
 * the number of distinct vertices in those cages.
 *
 * @param sys Frame handle.
 * @param cutoff Positive finite cutoff, Angstrom.
 * @param type Species filter.
 * @param graph ::SEAMS_GRAPH_CUTOFF or ::SEAMS_GRAPH_MUTUAL.
 * @param spec Signature text. Must be non-NULL.
 * @param n_cages Receives the cage count.
 * @param n_atoms Receives the vertex count.
 * @return ::SEAMS_OK, or ::SEAMS_ERR_NULL, ::SEAMS_ERR_STATE,
 *         ::SEAMS_ERR_RANGE, ::SEAMS_ERR_ARG, ::SEAMS_ERR_NOMEM,
 *         ::SEAMS_ERR_INTERNAL.
 * @par Handle access
 * Non-const only because the call clears the handle's error slot on
 * entry, so seams_handle_error() reports this call. It reads the frame
 * and keeps no scratch on the handle.
 * @par Ownership
 * The caller owns sys; the call borrows it. spec is read during the call
 * only. n_cages and n_atoms are caller storage that the call fills. The
 * library keeps no pointer to any argument.
 * @par Thread safety
 * Exclusive use of sys.
 * @par Example
 * @code
 * size_t n_cages = 0, n_atoms = 0;
 * seams_cages_signature(sys, 3.5, 0, SEAMS_GRAPH_CUTOFF, "512", &n_cages, &n_atoms);
 * @endcode
 */
SEAMS_API seams_status seams_cages_signature(seams_system *sys, double cutoff,
                                             int type, int graph, const char *spec,
                                             size_t *n_cages, size_t *n_atoms);

/**
 * @brief Per-atom ice shell label.
 *
 * labels[i] is 0 for ice Ih, 1 for ice Ic, 2 for stacking disorder, or
 * -1 when the shell does not match. A non-positive shell_rmsd selects
 * the library ceiling (the scale-free second-shell residual, 0.15).
 *
 * @param sys Frame handle.
 * @param shell_rmsd Positive scale-free residual ceiling, or a
 *        non-positive value for the library ceiling.
 * @param labels Caller buffer.
 * @param labels_len Must be at least the frame count.
 * @return ::SEAMS_OK, or ::SEAMS_ERR_NULL, ::SEAMS_ERR_STATE,
 *         ::SEAMS_ERR_RANGE, ::SEAMS_ERR_INTERNAL, ::SEAMS_ERR_NOMEM.
 * @par Handle access
 * Non-const only because the call clears the handle's error slot on
 * entry, so seams_handle_error() reports this call. It reads the frame
 * and keeps no scratch on the handle.
 * @par Ownership
 * The caller owns labels.
 * @par Thread safety
 * Exclusive use of sys.
 * @par Example
 * @code
 * seams_system_symmetry(sys, 0.0, labels, n);
 * @endcode
 */
SEAMS_API seams_status seams_system_symmetry(seams_system *sys, double shell_rmsd,
                                             int *labels, size_t labels_len);

/** @} */

/**
 * @defgroup capi_brine Ice/brine fronts
 * @brief Planar front, rejection coefficient, profile, and pockets.
 * @{
 */

/**
 * @brief Planar ice/brine front from an ice assignment.
 *
 * ice[i] is non-zero when atom i is ice. water_type 0 keeps every atom
 * whose type is not ion_type. ion_type 0 means the frame has no ions.
 * n_bins 0 selects 20 bins. axis is 0, 1, or 2. The ice fraction cut
 * is 0.5, and the bulk ice and liquid cuts are 0.85 and 0.15.
 *
 * @param sys Frame handle.
 * @param ice Per-atom ice flags, length n_ice.
 * @param n_ice Must equal the frame count.
 * @param water_type Species kept as water, or 0.
 * @param ion_type Species kept as ions, or 0 for none.
 * @param n_bins Bin count along axis, or 0 for 20.
 * @param axis 0, 1, or 2.
 * @param out Receives a new handle. On failure *out is NULL.
 * @return ::SEAMS_OK, or ::SEAMS_ERR_NULL, ::SEAMS_ERR_STATE,
 *         ::SEAMS_ERR_RANGE, ::SEAMS_ERR_ARG, ::SEAMS_ERR_NOMEM,
 *         ::SEAMS_ERR_INTERNAL.
 * @par Handle access
 * Non-const only because the call clears the handle's error slot on
 * entry, so seams_handle_error() reports this call. It reads the frame
 * and keeps no scratch on the handle.
 * @par Ownership
 * The caller frees *out with seams_brine_destroy().
 * @par Thread safety
 * Exclusive use of sys.
 * @par Example
 * @code
 * seams_brine *front = NULL;
 * if (seams_brine_analyse(sys, ice, n, 1, 2, 0, 2, &front) == SEAMS_OK) {
 *   seams_brine_destroy(front);
 * }
 * @endcode
 */
SEAMS_API seams_status seams_brine_analyse(seams_system *sys, const int *ice,
                                           size_t n_ice, int water_type, int ion_type,
                                           int n_bins, int axis, seams_brine **out);

/**
 * @brief Free a brine-front handle.
 * @param front Handle from seams_brine_analyse(), or NULL. NULL is ignored.
 * @par Ownership
 * Releases the handle and its profile.
 * @par Thread safety
 * The caller must not use front after this returns.
 * @par Example
 * @code
 * seams_brine_destroy(front);
 * @endcode
 */
SEAMS_API void seams_brine_destroy(seams_brine *front);

/**
 * @brief The two planar crossings of the ice slab.
 *
 * Lengths are Angstrom, in the frame's cell.
 *
 * @param front Brine handle.
 * @param z_enter Where the ice fraction crosses the cut on entering
 *        the ice along the positive axis, or NULL.
 * @param z_exit The crossing on leaving the ice, or NULL.
 * @param thickness Ice run along the axis, or NULL.
 * @param found Receives 1 when both crossings exist, otherwise 0, or NULL.
 * @return ::SEAMS_OK, or ::SEAMS_ERR_NULL when front is NULL or every
 *         output pointer is NULL.
 * @par Ownership
 * The caller owns front; the call borrows it. z_enter, z_exit, thickness and
 * found are caller storage that the call fills. The library keeps no pointer
 * to any argument.
 * @par Thread safety
 * Read-only: may run alongside other const calls on front, and not
 * alongside a non-const call on front.
 * @par Example
 * @code
 * double thick = 0;
 * int found = 0;
 * seams_brine_fronts(front, NULL, NULL, &thick, &found);
 * @endcode
 */
SEAMS_API seams_status seams_brine_fronts(const seams_brine *front, double *z_enter,
                                          double *z_exit, double *thickness, int *found);

/**
 * @brief Molal rejection coefficient k = c_ice / c_liquid.
 *
 * Molality is formula units per kilogram of water.
 *
 * @param front Brine handle.
 * @param c_ice Bulk ice molality, or NULL.
 * @param c_liquid Bulk liquid molality, or NULL.
 * @param k The ratio, or NULL. Unchanged when the ratio is undefined.
 * @param defined Receives 1 when the liquid bulk contains salt or the
 *        ice bulk does not, otherwise 0, or NULL.
 * @return ::SEAMS_OK, or ::SEAMS_ERR_NULL when front is NULL or every
 *         output pointer is NULL.
 * @par Ownership
 * The caller owns front; the call borrows it. c_ice, c_liquid, k and defined
 * are caller storage that the call fills. The library keeps no pointer to any
 * argument.
 * @par Thread safety
 * Read-only: may run alongside other const calls on front, and not
 * alongside a non-const call on front.
 * @par Example
 * @code
 * double k = 0;
 * int defined = 0;
 * seams_brine_rejection(front, NULL, NULL, &k, &defined);
 * @endcode
 */
SEAMS_API seams_status seams_brine_rejection(const seams_brine *front, double *c_ice,
                                             double *c_liquid, double *k, int *defined);

/**
 * @brief Concentration profile along the front axis.
 *
 * z_mid is Angstrom. molality is mol/kg. ice_fraction is dimensionless.
 *
 * @param front Brine handle.
 * @param z_mid Bin centres, or NULL to skip the column.
 * @param molality Bin molalities, or NULL to skip the column.
 * @param ice_fraction Bin ice fractions, or NULL to skip the column.
 * @param cap Capacity of each non-NULL column.
 * @param n_out Receives the bin count.
 * @return ::SEAMS_OK when every non-NULL column was filled, or when every
 *         column pointer was NULL and n_out was written.
 *         ::SEAMS_ERR_RANGE when a non-NULL column is short.
 *         ::SEAMS_ERR_NULL when front or n_out is NULL.
 * @par Ownership
 * The caller owns the columns.
 * @par Thread safety
 * Read-only: may run alongside other const calls on front, and not
 * alongside a non-const call on front.
 * @par Example
 * @code
 * size_t bins = 0;
 * seams_brine_profile(front, NULL, NULL, NULL, 0, &bins);
 * @endcode
 */
SEAMS_API seams_status seams_brine_profile(const seams_brine *front, double *z_mid,
                                           double *molality, double *ice_fraction,
                                           size_t cap, size_t *n_out);

/**
 * @brief Brine pockets trapped in the ice, excluding the liquid reservoir.
 * @param front Brine handle.
 * @param n_pockets Receives the pocket count, or NULL.
 * @param reservoir_waters Receives the water count of the liquid
 *        reservoir, or NULL.
 * @return ::SEAMS_OK, or ::SEAMS_ERR_NULL when front is NULL or both
 *         output pointers are NULL.
 * @par Ownership
 * The caller owns front; the call borrows it. n_pockets and reservoir_waters
 * are caller storage that the call fills. The library keeps no pointer to any
 * argument.
 * @par Thread safety
 * Read-only: may run alongside other const calls on front, and not
 * alongside a non-const call on front.
 * @par Example
 * @code
 * int pockets = 0;
 * seams_brine_pockets(front, &pockets, NULL);
 * @endcode
 */
SEAMS_API seams_status seams_brine_pockets(const seams_brine *front, int *n_pockets,
                                           int *reservoir_waters);

/** @} */

/**
 * @defgroup capi_more Further analyses
 * @brief Aliases, k-nearest lists, RDF, hydrogen bonds, density, fingerprints, sites.
 * @{
 */

/**
 * @brief Hexagonal and double-diamond cages. Same search as seams_cages_find().
 * @param sys Frame handle.
 * @param cutoff Positive finite cutoff, Angstrom.
 * @param type Species filter.
 * @param graph ::SEAMS_GRAPH_CUTOFF, ::SEAMS_GRAPH_MUTUAL, or ::SEAMS_GRAPH_SEEDED.
 * @param k Neighbour count for the seeded graph. Zero selects 4.
 * @param complete Non-zero runs ring-adjacent completion on seeded labels.
 * @param out Receives a ::seams_cage_set. On failure *out is NULL.
 * @return Same codes as seams_cages_find().
 * @par Handle access
 * Non-const only because the call clears the handle's error slot on
 * entry, so seams_handle_error() reports this call. It reads the frame
 * and keeps no scratch on the handle.
 * @par Ownership
 * The caller frees *out with seams_cages_destroy().
 * @par Thread safety
 * Exclusive use of sys.
 * @par Example
 * @code
 * seams_cage_set *cages = NULL;
 * seams_cages(sys, 3.5, 0, SEAMS_GRAPH_SEEDED, 4, 1, &cages);
 * seams_cages_destroy(cages);
 * @endcode
 */
SEAMS_API seams_status seams_cages(seams_system *sys, double cutoff, int type, int graph,
                                   int k, int complete, seams_cage_set **out);

/**
 * @brief Closed cages of one ring census. Same search as seams_cages_signature().
 * @param sys Frame handle.
 * @param cutoff Positive finite cutoff, Angstrom.
 * @param type Species filter.
 * @param graph ::SEAMS_GRAPH_CUTOFF or ::SEAMS_GRAPH_MUTUAL.
 * @param spec Signature text such as "512" or "hc". Must be non-NULL.
 * @param n_cages Receives the cage count.
 * @param n_atoms Receives the vertex count.
 * @return Same codes as seams_cages_signature().
 * @par Handle access
 * Non-const only because the call clears the handle's error slot on
 * entry, so seams_handle_error() reports this call. It reads the frame
 * and keeps no scratch on the handle.
 * @par Ownership
 * The caller owns sys; the call borrows it. spec is read during the call
 * only. n_cages and n_atoms are caller storage that the call fills. The
 * library keeps no pointer to any argument.
 * @par Thread safety
 * Exclusive use of sys.
 * @par Example
 * @code
 * size_t n_cages = 0, n_atoms = 0;
 * seams_signature_cages(sys, 3.5, 0, SEAMS_GRAPH_CUTOFF, "hc", &n_cages, &n_atoms);
 * @endcode
 */
SEAMS_API seams_status seams_signature_cages(seams_system *sys, double cutoff, int type,
                                             int graph, const char *spec, size_t *n_cages,
                                             size_t *n_atoms);

/**
 * @brief Local Steinhardt q_l on the cutoff shell, and the average on request.
 *
 * This is seams_ql() with the cutoff graph.
 *
 * @param sys Frame handle.
 * @param cutoff Positive finite cutoff, Angstrom.
 * @param type Species filter.
 * @param order_l 3, 4, 6, 8, or 12.
 * @param ql Local q_l, or NULL.
 * @param ql_len Length of ql. Ignored when ql is NULL.
 * @param qlbar Averaged q_l, or NULL to skip the Lechner-Dellago average.
 * @param qlbar_len Length of qlbar. Ignored when qlbar is NULL.
 * @return Same codes as seams_ql().
 * @par Handle access
 * Non-const only because the call clears the handle's error slot on
 * entry; it reads the frame and nothing else on the handle. The q_lm work
 * arrays belong to the calling thread, so calls on different handles may
 * run at the same time.
 * @par Ownership
 * The caller owns both arrays.
 * @par Thread safety
 * Exclusive use of sys. Any q call on another handle may run beside it.
 * @par Example
 * @code
 * seams_steinhardt(sys, 3.5, 0, 6, q6, n, NULL, 0);
 * @endcode
 */
SEAMS_API seams_status seams_steinhardt(seams_system *sys, double cutoff, int type,
                                        int order_l, double *ql, size_t ql_len,
                                        double *qlbar, size_t qlbar_len);

/**
 * @brief k-nearest neighbour list.
 * @param sys Frame handle.
 * @param k Positive neighbour count.
 * @param cutoff Nominal positive cutoff, Angstrom. Candidates are taken
 *        inside cutoff + 1.5.
 * @param type Species filter. Zero and negative keep every atom.
 * @param mutual Non-zero keeps the mutual graph. Zero keeps the one-way
 *        k-nearest graph.
 * @param out Receives a neighbour handle. On failure *out is NULL.
 * @return ::SEAMS_OK, or ::SEAMS_ERR_NULL, ::SEAMS_ERR_STATE,
 *         ::SEAMS_ERR_RANGE, ::SEAMS_ERR_ARG, ::SEAMS_ERR_NOMEM,
 *         ::SEAMS_ERR_INTERNAL.
 * @par Handle access
 * Non-const only because the call clears the handle's error slot on
 * entry, so seams_handle_error() reports this call. It reads the frame
 * and keeps no scratch on the handle.
 * @par Ownership
 * The caller frees *out with seams_neighbours_destroy().
 * @par Thread safety
 * Exclusive use of sys.
 * @par Example
 * @code
 * seams_neighbours *nb = NULL;
 * seams_neighbours_knn(sys, 4, 3.5, 0, 1, &nb);
 * seams_neighbours_destroy(nb);
 * @endcode
 */
SEAMS_API seams_status seams_neighbours_knn(seams_system *sys, int k, double cutoff,
                                            int type, int mutual, seams_neighbours **out);

/**
 * @brief Partial radial distribution.
 *
 * r is Angstrom. g is dimensionless. volume is cubic Angstrom.
 *
 * @param sys Frame handle.
 * @param rmax Positive finite outer radius, Angstrom.
 * @param bins Bin count. Zero selects about rmax / 0.1.
 * @param type_i First species. Zero and negative keep every atom of a
 *        single-species frame, or every atom as one species otherwise.
 * @param type_j Second species. Zero and negative reuse the type_i species.
 * @param r Bin centres, or NULL to query.
 * @param g g(r) per bin, or NULL to query.
 * @param count Pair count per bin, or NULL to query.
 * @param cap Capacity of each non-NULL column.
 * @param n_out Receives the bin count.
 * @param volume Receives the cell volume, or NULL.
 * @return ::SEAMS_OK, or ::SEAMS_ERR_NULL, ::SEAMS_ERR_STATE,
 *         ::SEAMS_ERR_RANGE, ::SEAMS_ERR_NOMEM, ::SEAMS_ERR_INTERNAL.
 * @par Handle access
 * Non-const only because the call clears the handle's error slot on
 * entry, so seams_handle_error() reports this call. It reads the frame
 * and keeps no scratch on the handle.
 * @par Ownership
 * The caller owns the columns.
 * @par Thread safety
 * Exclusive use of sys.
 * @par Example
 * @code
 * size_t bins = 0;
 * seams_rdf(sys, 8.0, 80, 0, 0, NULL, NULL, NULL, 0, &bins, NULL);
 * @endcode
 */
SEAMS_API seams_status seams_rdf(seams_system *sys, double rmax, int bins, int type_i,
                                 int type_j, double *r, double *g, int *count, size_t cap,
                                 size_t *n_out, double *volume);

/**
 * @brief Running coordination number along the partial RDF.
 * @param sys Frame handle.
 * @param rmax Positive finite outer radius, Angstrom.
 * @param bins Bin count. Zero selects about rmax / 0.1.
 * @param type_i Central species, as in seams_rdf().
 * @param type_j Counted species, as in seams_rdf().
 * @param r Bin centres, Angstrom, or NULL to query.
 * @param cn Running coordination number, or NULL to query.
 * @param cap Capacity of each non-NULL column.
 * @param n_out Receives the bin count.
 * @return ::SEAMS_OK, or ::SEAMS_ERR_NULL, ::SEAMS_ERR_STATE,
 *         ::SEAMS_ERR_RANGE, ::SEAMS_ERR_NOMEM, ::SEAMS_ERR_INTERNAL.
 * @par Handle access
 * Non-const only because the call clears the handle's error slot on
 * entry, so seams_handle_error() reports this call. It reads the frame
 * and keeps no scratch on the handle.
 * @par Ownership
 * The caller owns the columns.
 * @par Thread safety
 * Exclusive use of sys.
 * @par Example
 * @code
 * double r[80], cn[80];
 * size_t bins = 0;
 * seams_running_cn(sys, 8.0, 80, 0, 0, r, cn, 80, &bins);
 * @endcode
 */
SEAMS_API seams_status seams_running_cn(seams_system *sys, double rmax, int bins, int type_i,
                                        int type_j, double *r, double *cn, size_t cap,
                                        size_t *n_out);

/**
 * @brief Coordination number integrated to rmax.
 * @param sys Frame handle.
 * @param rmax Positive finite outer radius, Angstrom.
 * @param bins Bin count. Zero selects about rmax / 0.1.
 * @param type_i Central species, as in seams_rdf().
 * @param type_j Counted species, as in seams_rdf().
 * @param cn Receives the coordination number.
 * @return ::SEAMS_OK, or ::SEAMS_ERR_NULL, ::SEAMS_ERR_STATE,
 *         ::SEAMS_ERR_RANGE, ::SEAMS_ERR_NOMEM, ::SEAMS_ERR_INTERNAL.
 * @par Handle access
 * Non-const only because the call clears the handle's error slot on
 * entry, so seams_handle_error() reports this call. It reads the frame
 * and keeps no scratch on the handle.
 * @par Ownership
 * The caller owns sys; the call borrows it. cn is caller storage that the
 * call fills. The library keeps no pointer to any argument.
 * @par Thread safety
 * Exclusive use of sys.
 * @par Example
 * @code
 * double cn = 0;
 * seams_coordination(sys, 3.5, 0, 1, 1, &cn);
 * @endcode
 */
SEAMS_API seams_status seams_coordination(seams_system *sys, double rmax, int bins,
                                          int type_i, int type_j, double *cn);

/**
 * @brief Hydrogen-bond count on one frame.
 * @param sys Frame handle.
 * @param cutoff Positive finite donor-acceptor cutoff, Angstrom.
 * @param oxygen_type Heavy-atom species. Zero keeps every atom that is
 *        not hydrogen_type.
 * @param hydrogen_type Positive hydrogen species.
 * @param dist_cutoff Acceptor-H distance, Angstrom.
 * @param angle_cutoff O-O-H angle, degrees.
 * @param n_bonds Receives the number of hydrogen bonds.
 * @return ::SEAMS_OK, or ::SEAMS_ERR_NULL, ::SEAMS_ERR_STATE,
 *         ::SEAMS_ERR_RANGE, ::SEAMS_ERR_ARG, ::SEAMS_ERR_NOMEM,
 *         ::SEAMS_ERR_INTERNAL.
 * @par Handle access
 * Non-const only because the call clears the handle's error slot on
 * entry, so seams_handle_error() reports this call. It reads the frame
 * and keeps no scratch on the handle.
 * @par Ownership
 * The caller owns sys; the call borrows it. n_bonds is caller storage that
 * the call fills. The library keeps no pointer to any argument.
 * @par Thread safety
 * Exclusive use of sys.
 * @par Example
 * @code
 * int n_bonds = 0;
 * seams_hbonds(sys, 3.5, 1, 2, 2.4, 30.0, &n_bonds);
 * @endcode
 */
SEAMS_API seams_status seams_hbonds(seams_system *sys, double cutoff, int oxygen_type,
                                    int hydrogen_type, double dist_cutoff,
                                    double angle_cutoff, int *n_bonds);

/**
 * @brief Number density along one cell axis.
 *
 * z is Angstrom. rho is number density, atoms per cubic Angstrom.
 *
 * @param sys Frame handle.
 * @param type Species filter. Zero keeps every atom.
 * @param bins Bin count. Zero selects about one bin per 0.1 Angstrom.
 * @param axis 0, 1, or 2.
 * @param z Bin centres, or NULL to query.
 * @param rho Density per bin, or NULL to query.
 * @param cap Capacity of each non-NULL column.
 * @param n_out Receives the bin count.
 * @return ::SEAMS_OK, or ::SEAMS_ERR_NULL, ::SEAMS_ERR_STATE,
 *         ::SEAMS_ERR_RANGE, ::SEAMS_ERR_NOMEM, ::SEAMS_ERR_INTERNAL.
 * @par Handle access
 * Non-const only because the call clears the handle's error slot on
 * entry, so seams_handle_error() reports this call. It reads the frame
 * and keeps no scratch on the handle.
 * @par Ownership
 * The caller owns the columns.
 * @par Thread safety
 * Exclusive use of sys.
 * @par Example
 * @code
 * double z[50], rho[50];
 * size_t bins = 0;
 * seams_density_axis(sys, 0, 50, 2, z, rho, 50, &bins);
 * @endcode
 */
SEAMS_API seams_status seams_density_axis(seams_system *sys, int type, int bins, int axis,
                                          double *z, double *rho, size_t cap, size_t *n_out);

/**
 * @brief Topological fingerprint of the bonded graph.
 * @param sys Frame handle.
 * @param cutoff Positive finite cutoff, Angstrom.
 * @param type Species filter.
 * @param k Neighbour count for the mutual graph. Zero selects 4.
 * @param graph ::SEAMS_GRAPH_CUTOFF or ::SEAMS_GRAPH_MUTUAL.
 * @param hops Neighbourhood depth. Zero selects 2.
 * @param colour_types Non-zero colours vertices by LAMMPS type.
 * @param out Receives a ::seams_topo. On failure *out is NULL.
 * @return ::SEAMS_OK, or ::SEAMS_ERR_NULL, ::SEAMS_ERR_STATE,
 *         ::SEAMS_ERR_RANGE, ::SEAMS_ERR_ARG, ::SEAMS_ERR_NOMEM,
 *         ::SEAMS_ERR_INTERNAL.
 * @par Handle access
 * Non-const only because the call clears the handle's error slot on
 * entry, so seams_handle_error() reports this call. It reads the frame
 * and keeps no scratch on the handle.
 * @par Ownership
 * The caller frees *out with seams_fingerprint_destroy().
 * @par Thread safety
 * Exclusive use of sys.
 * @par Example
 * @code
 * seams_topo *fp = NULL;
 * seams_fingerprint(sys, 3.5, 0, 4, SEAMS_GRAPH_CUTOFF, 2, 0, &fp);
 * seams_fingerprint_destroy(fp);
 * @endcode
 */
SEAMS_API seams_status seams_fingerprint(seams_system *sys, double cutoff, int type, int k,
                                         int graph, int hops, int colour_types,
                                         seams_topo **out);

/**
 * @brief Free a fingerprint.
 * @param fp Handle from seams_fingerprint(), or NULL. NULL is ignored.
 * @par Ownership
 * Releases the handle and its key strings.
 * @par Thread safety
 * The caller must not use fp after this returns.
 * @par Example
 * @code
 * seams_fingerprint_destroy(fp);
 * @endcode
 */
SEAMS_API void seams_fingerprint_destroy(seams_topo *fp);

/**
 * @brief Frame key text.
 * @param fp Fingerprint handle.
 * @param buf Caller buffer of cap bytes, or NULL to query.
 * @param cap Capacity of buf, including the terminating NUL.
 * @param n_out Receives the key length without the terminating NUL.
 * @return ::SEAMS_OK, or ::SEAMS_ERR_NULL, ::SEAMS_ERR_RANGE when buf is
 *         short.
 * @par Ownership
 * The caller owns buf. The copy is NUL-terminated.
 * @par Thread safety
 * Read-only: may run alongside other const calls on fp, and not
 * alongside a non-const call on fp.
 * @par Example
 * @code
 * size_t len = 0;
 * seams_fingerprint_key(fp, NULL, 0, &len);
 * char *key = malloc(len + 1);
 * seams_fingerprint_key(fp, key, len + 1, &len);
 * @endcode
 */
SEAMS_API seams_status seams_fingerprint_key(const seams_topo *fp, char *buf, size_t cap,
                                             size_t *n_out);

/**
 * @brief Local key of one atom. Same buffer rule as seams_fingerprint_key().
 * @param fp Fingerprint handle.
 * @param index Atom index in the frame.
 * @param buf Caller buffer of cap bytes, or NULL to query.
 * @param cap Capacity of buf, including the terminating NUL.
 * @param n_out Receives the key length without the terminating NUL.
 * @return ::SEAMS_OK, or ::SEAMS_ERR_NULL, ::SEAMS_ERR_RANGE for a bad
 *         index or a short buffer.
 * @par Ownership
 * The caller owns buf.
 * @par Thread safety
 * Read-only: may run alongside other const calls on fp, and not
 * alongside a non-const call on fp.
 * @par Example
 * @code
 * char key[256];
 * size_t len = 0;
 * seams_fingerprint_atom(fp, 0, key, sizeof key, &len);
 * @endcode
 */
SEAMS_API seams_status seams_fingerprint_atom(const seams_topo *fp, int index, char *buf,
                                              size_t cap, size_t *n_out);

/**
 * @brief Primitive-ring census stored on the fingerprint.
 * @param fp Fingerprint handle.
 * @param counts Caller buffer. counts[k] receives the number of rings of
 *        size k. Sizes past the buffer are dropped.
 * @param counts_len Length of counts. Must be positive.
 * @return ::SEAMS_OK, or ::SEAMS_ERR_NULL, ::SEAMS_ERR_RANGE when
 *         counts_len is 0.
 * @par Ownership
 * The caller owns counts.
 * @par Thread safety
 * Read-only: may run alongside other const calls on fp, and not
 * alongside a non-const call on fp.
 * @par Example
 * @code
 * int census[13];
 * seams_fingerprint_rings(fp, census, 13);
 * @endcode
 */
SEAMS_API seams_status seams_fingerprint_rings(const seams_topo *fp, int *counts,
                                               size_t counts_len);

/**
 * @brief Mutual nearest unlike contacts of a site table.
 *
 * Cation and anion counts are molecules in the ion cloud.
 *
 * @param sys Frame handle.
 * @param site_spec Text such as "1=cationHead,2=anion". Must be non-NULL.
 * @param n_pairs Receives the contact count, or NULL.
 * @param n_cation Receives the cation molecule count, or NULL.
 * @param n_anion Receives the anion molecule count, or NULL. At least
 *        one of the three outputs must be non-NULL.
 * @return ::SEAMS_OK, or ::SEAMS_ERR_NULL, ::SEAMS_ERR_STATE,
 *         ::SEAMS_ERR_ARG, ::SEAMS_ERR_NOMEM, ::SEAMS_ERR_INTERNAL.
 * @par Handle access
 * Non-const only because the call clears the handle's error slot on
 * entry, so seams_handle_error() reports this call. It reads the frame
 * and keeps no scratch on the handle.
 * @par Ownership
 * The caller owns sys; the call borrows it. site_spec is read during the call
 * only. n_pairs, n_cation and n_anion are caller storage that the call fills.
 * The library keeps no pointer to any argument.
 * @par Thread safety
 * Exclusive use of sys.
 * @par Example
 * @code
 * size_t pairs = 0;
 * seams_contact_pairs(sys, "1=cationHead,2=anion", &pairs, NULL, NULL);
 * @endcode
 */
SEAMS_API seams_status seams_contact_pairs(seams_system *sys, const char *site_spec,
                                           size_t *n_pairs, int *n_cation, int *n_anion);

/**
 * @brief Largest domain of a polar or apolar site subset.
 * @param sys Frame handle.
 * @param cutoff Positive finite cluster cutoff, Angstrom.
 * @param site_spec Site table text. Must be non-NULL.
 * @param subset 0 is apolar, 1 is polar.
 * @param n_subset Receives the number of sites in the subset, or NULL.
 * @param largest Receives the size of the largest domain, or NULL.
 * @param percolation Receives the largest domain size over the subset
 *        size, or NULL. At least one output must be non-NULL.
 * @return ::SEAMS_OK, or ::SEAMS_ERR_NULL, ::SEAMS_ERR_STATE,
 *         ::SEAMS_ERR_RANGE, ::SEAMS_ERR_ARG, ::SEAMS_ERR_NOMEM,
 *         ::SEAMS_ERR_INTERNAL.
 * @par Handle access
 * Non-const only because the call clears the handle's error slot on
 * entry, so seams_handle_error() reports this call. It reads the frame
 * and keeps no scratch on the handle.
 * @par Ownership
 * The caller owns sys; the call borrows it. site_spec is read during the call
 * only. n_subset, largest and percolation are caller storage that the call
 * fills. The library keeps no pointer to any argument.
 * @par Thread safety
 * Exclusive use of sys.
 * @par Example
 * @code
 * int largest = 0;
 * seams_domains(sys, 5.0, "1=cationHead,2=anion,3=tail", 0, NULL, &largest, NULL);
 * @endcode
 */
SEAMS_API seams_status seams_domains(seams_system *sys, double cutoff, const char *site_spec,
                                     int subset, int *n_subset, int *largest,
                                     double *percolation);

/**
 * @brief Perfect and imperfect prism counts by ring size.
 *
 * Rings come from the mutual four-nearest graph inside cutoff + 1.5.
 *
 * @param sys Frame handle.
 * @param cutoff Positive finite cutoff, Angstrom.
 * @param type Species filter.
 * @param perfect Caller buffer. perfect[k] receives the perfect k-prisms.
 * @param imperfect Caller buffer. imperfect[k] receives the deformed
 *        k-prisms.
 * @param len Length of each buffer. At least 7, so that sizes 3 through
 *        6 are filled.
 * @return ::SEAMS_OK, or ::SEAMS_ERR_NULL, ::SEAMS_ERR_STATE,
 *         ::SEAMS_ERR_RANGE, ::SEAMS_ERR_NOMEM, ::SEAMS_ERR_INTERNAL.
 * @par Handle access
 * Non-const only because the call clears the handle's error slot on
 * entry, so seams_handle_error() reports this call. It reads the frame
 * and keeps no scratch on the handle.
 * @par Ownership
 * The caller owns both buffers.
 * @par Thread safety
 * Exclusive use of sys.
 * @par Example
 * @code
 * int perfect[7], imperfect[7];
 * seams_prisms(sys, 3.5, 0, perfect, imperfect, 7);
 * @endcode
 */
SEAMS_API seams_status seams_prisms(seams_system *sys, double cutoff, int type, int *perfect,
                                    int *imperfect, size_t len);

/**
 * @brief Primitive-ring histogram of a monolayer.
 * @param sys Frame handle.
 * @param cutoff Positive finite cutoff, Angstrom.
 * @param type Species filter.
 * @param max_depth Largest ring size, at least 3.
 * @param sheet_area Sheet area, square Angstrom. Must be positive.
 * @param counts counts[k] is the number of rings of size k.
 * @param counts_len Must be at least max_depth + 1.
 * @return ::SEAMS_OK, or ::SEAMS_ERR_NULL, ::SEAMS_ERR_STATE,
 *         ::SEAMS_ERR_RANGE, ::SEAMS_ERR_NOMEM, ::SEAMS_ERR_INTERNAL.
 * @par Handle access
 * Non-const only because the call clears the handle's error slot on
 * entry, so seams_handle_error() reports this call. It reads the frame
 * and keeps no scratch on the handle.
 * @par Ownership
 * The caller owns counts.
 * @par Thread safety
 * Exclusive use of sys.
 * @par Example
 * @code
 * int counts[7];
 * seams_monolayer(sys, 3.5, 0, 6, 900.0, counts, 7);
 * @endcode
 */
SEAMS_API seams_status seams_monolayer(seams_system *sys, double cutoff, int type,
                                       int max_depth, double sheet_area, int *counts,
                                       size_t counts_len);

/**
 * @brief Per-atom site kind from a site-table spec.
 * @param sys Frame handle.
 * @param spec Site table text such as "1=cationHead,2=anion". Must be
 *        non-NULL.
 * @param kinds One ::seams_site_kind per atom.
 * @param len Length of kinds. Must be at least the frame count.
 * @return ::SEAMS_OK, or ::SEAMS_ERR_NULL, ::SEAMS_ERR_STATE,
 *         ::SEAMS_ERR_RANGE, ::SEAMS_ERR_ARG, ::SEAMS_ERR_NOMEM,
 *         ::SEAMS_ERR_INTERNAL.
 * @par Handle access
 * Non-const only because the call clears the handle's error slot on
 * entry, so seams_handle_error() reports this call. It reads the frame
 * and keeps no scratch on the handle.
 * @par Ownership
 * The caller owns kinds.
 * @par Thread safety
 * Exclusive use of sys.
 * @par Example
 * @code
 * seams_site_table(sys, "1=cationHead,2=anion", kinds, n);
 * @endcode
 */
SEAMS_API seams_status seams_site_table(seams_system *sys, const char *spec, int *kinds,
                                        size_t len);

/** @} */

/**
 * @defgroup capi_sized Status-returning and length-checked entry points
 * @brief Preferred spellings for new callers. Each one reports errors
 *        through ::seams_status and takes a size_t length for every
 *        caller buffer.
 * @{
 */

/**
 * @brief Thread-local detail of the last failing call on this thread.
 *
 * The string is empty after a successful call on this thread. Prefer
 * seams_handle_error() when the failing call had a handle: that copy
 * survives a later call on another handle. seams_last_error() is the
 * ABI 1 spelling and returns the same string.
 *
 * @return Pointer to a thread-local string. Never NULL.
 * @par Ownership
 * Owned by the library. Valid until the next call on this thread that
 * sets or clears the error.
 * @par Thread safety
 * The storage is thread-local.
 * @par Example
 * @code
 * seams_status rc = seams_system_natoms(NULL, NULL);
 * fprintf(stderr, "%s: %s\n", seams_status_string(rc), seams_last_error_message());
 * @endcode
 */
SEAMS_API const char *seams_last_error_message(void);

/**
 * @brief Number of points in the frame.
 * @param sys Frame handle.
 * @param n Receives the atom count.
 * @return ::SEAMS_OK, ::SEAMS_ERR_NULL when sys or n is NULL, or
 *         ::SEAMS_ERR_STATE when sys holds no frame. *n is unchanged on
 *         failure.
 * @par Ownership
 * The caller owns sys; the call borrows it. n is caller storage that the call
 * fills. The library keeps no pointer to any argument.
 * @par Thread safety
 * Read-only: may run alongside other const calls on sys, and not
 * alongside a non-const call on sys.
 * @par Example
 * @code
 * size_t n = 0;
 * if (seams_system_natoms(sys, &n) == SEAMS_OK) { labels = malloc(n * sizeof *labels); }
 * @endcode
 */
SEAMS_API seams_status seams_system_natoms(const seams_system *sys, size_t *n);

/**
 * @brief Number of neighbour rows.
 * @param nb Neighbour handle.
 * @param n Receives the row count.
 * @return ::SEAMS_OK, or ::SEAMS_ERR_NULL when nb or n is NULL.
 * @par Ownership
 * The caller owns nb; the call borrows it. n is caller storage that the call
 * fills. The library keeps no pointer to any argument.
 * @par Thread safety
 * Read-only: may run alongside other const calls on nb, and not
 * alongside a non-const call on nb.
 * @par Example
 * @code
 * size_t rows = 0;
 * seams_neighbours_nrows(nb, &rows);
 * @endcode
 */
SEAMS_API seams_status seams_neighbours_nrows(const seams_neighbours *nb, size_t *n);

/**
 * @brief Number of primitive rings.
 * @param rings Ring handle.
 * @param n Receives the ring count.
 * @return ::SEAMS_OK, or ::SEAMS_ERR_NULL when rings or n is NULL.
 * @par Ownership
 * The caller owns rings; the call borrows it. n is caller storage that the
 * call fills. The library keeps no pointer to any argument.
 * @par Thread safety
 * Read-only: may run alongside other const calls on rings, and not
 * alongside a non-const call on rings.
 * @par Example
 * @code
 * size_t nr = 0;
 * seams_rings_nrings(rings, &nr);
 * @endcode
 */
SEAMS_API seams_status seams_rings_nrings(const seams_rings *rings, size_t *n);

/**
 * @brief CHILL+ labels into a buffer of known length.
 *
 * rule selects the q_lm shell. ::SEAMS_CHILL_PUBLISHED, the default, is
 * Nguyen and Molinero equation 1: the four closest neighbours.
 * ::SEAMS_CHILL_OVITO is every neighbour inside the cutoff, which is the
 * OVITO Chill+ modifier. Labels match ::seams_label.
 *
 * @param sys Frame handle.
 * @param cutoff Positive finite cutoff, Angstrom.
 * @param type Species filter. Zero and negative keep every atom.
 * @param graph ::SEAMS_GRAPH_CUTOFF or ::SEAMS_GRAPH_MUTUAL.
 * @param rule ::SEAMS_CHILL_PUBLISHED or ::SEAMS_CHILL_OVITO.
 * @param labels Caller buffer of labels_len ints.
 * @param labels_len Must be at least the frame count.
 * @param counts9 Nine ints that receive the population of each label,
 *        or NULL.
 * @return ::SEAMS_OK; ::SEAMS_ERR_NULL for a NULL sys or labels;
 *         ::SEAMS_ERR_STATE when sys holds no frame; ::SEAMS_ERR_RANGE
 *         for a short labels buffer or a non-positive or non-finite
 *         cutoff; ::SEAMS_ERR_ARG for an unknown rule or graph;
 *         ::SEAMS_ERR_NOMEM; ::SEAMS_ERR_INTERNAL.
 * @par Handle access
 * Non-const: writes each atom's CHILL+ label into the stored frame and
 * clears the handle's error slot on entry. Coordinates, cell, types, and
 * ids are not changed.
 * @par Ownership
 * The caller owns labels and counts9.
 * @par Thread safety
 * Exclusive use of sys.
 * @par Example
 * @code
 * int counts[9];
 * seams_chill_plus_n(sys, 3.5, 0, SEAMS_GRAPH_CUTOFF, SEAMS_CHILL_PUBLISHED,
 *                    labels, n, counts);
 * @endcode
 */
SEAMS_API seams_status seams_chill_plus_n(seams_system *sys, double cutoff, int type,
                                          int graph, int rule, int *labels,
                                          size_t labels_len, int counts9[9]);

/**
 * @brief Local q4 and q6 on the cutoff shell into buffers of known length.
 * @param sys Frame handle.
 * @param cutoff Positive finite cutoff, Angstrom.
 * @param type Species filter.
 * @param q4 Local q4, or NULL.
 * @param q4_len Length of q4. Ignored when q4 is NULL.
 * @param q6 Local q6, or NULL. q4 and q6 may not both be NULL.
 * @param q6_len Length of q6. Ignored when q6 is NULL.
 * @return ::SEAMS_OK, or ::SEAMS_ERR_NULL, ::SEAMS_ERR_STATE,
 *         ::SEAMS_ERR_RANGE for a short buffer or a bad cutoff,
 *         ::SEAMS_ERR_NOMEM, ::SEAMS_ERR_INTERNAL.
 * @par Handle access
 * Non-const only because the call clears the handle's error slot on
 * entry; it reads the frame and nothing else on the handle. The q_lm work
 * arrays belong to the calling thread, so calls on different handles may
 * run at the same time.
 * @par Ownership
 * The caller owns q4 and q6.
 * @par Thread safety
 * Exclusive use of sys. Any q call on another handle may run beside it.
 * @par Example
 * @code
 * seams_q4q6_n(sys, 3.5, 0, q4, n, q6, n);
 * @endcode
 */
SEAMS_API seams_status seams_q4q6_n(seams_system *sys, double cutoff, int type, double *q4,
                                    size_t q4_len, double *q6, size_t q6_len);

/**
 * @brief Size and symmetry of each ring into buffers of known length.
 * @param rings Ring handle.
 * @param sizes One int per ring, or NULL.
 * @param sizes_len Length of sizes. Ignored when sizes is NULL.
 * @param symmetry One ::seams_ring_symmetry code per ring, or NULL. The
 *        two pointers may not both be NULL.
 * @param symmetry_len Length of symmetry. Ignored when symmetry is NULL.
 * @return ::SEAMS_OK, or ::SEAMS_ERR_NULL, ::SEAMS_ERR_RANGE when a
 *         buffer is shorter than the ring count, ::SEAMS_ERR_NOMEM,
 *         ::SEAMS_ERR_INTERNAL.
 * @par Ownership
 * The caller owns both arrays.
 * @par Thread safety
 * Read-only: may run alongside other const calls on rings, and not
 * alongside a non-const call on rings.
 * @par Example
 * @code
 * size_t nr = 0;
 * seams_rings_nrings(rings, &nr);
 * int *sym = malloc(nr * sizeof *sym);
 * seams_rings_labels_n(rings, NULL, 0, sym, nr);
 * @endcode
 */
SEAMS_API seams_status seams_rings_labels_n(const seams_rings *rings, int *sizes,
                                            size_t sizes_len, int *symmetry,
                                            size_t symmetry_len);

/**
 * @brief Count primitive rings by size into a buffer of known length.
 * @param rings Ring handle.
 * @param counts Caller buffer. counts[k] receives the number of rings
 *        of size k. Sizes outside the buffer are ignored.
 * @param counts_len Length of counts. Must be positive.
 * @return ::SEAMS_OK, or ::SEAMS_ERR_NULL, ::SEAMS_ERR_RANGE when
 *         counts_len is 0.
 * @par Ownership
 * The caller owns counts.
 * @par Thread safety
 * Read-only: may run alongside other const calls on rings, and not
 * alongside a non-const call on rings.
 * @par Example
 * @code
 * int hist[13];
 * seams_rings_histogram_n(rings, hist, 13);
 * @endcode
 */
SEAMS_API seams_status seams_rings_histogram_n(const seams_rings *rings, int *counts,
                                               size_t counts_len);

/** @} */

/**
 * @defgroup capi_frame_io Frame reload and read-back
 * @brief Overwrite coordinates in place and copy the stored frame out.
 *
 * There is no borrowed-pointer load. Every analysis reads the library's
 * per-atom records, so a view onto caller memory would be copied on each
 * analysis call instead of once per frame. A trajectory loop calls
 * seams_system_set_atoms() once and seams_system_set_positions() for each
 * later frame, which writes only the coordinates.
 * @{
 */

/**
 * @brief Overwrite the coordinates, and optionally the cell, of the stored frame.
 *
 * Types and ids stay as stored. Per-atom labels written by an earlier
 * seams_chill() or seams_chill_plus() are reset. On failure the frame is
 * unchanged.
 *
 * By default the call checks every coordinate first and then writes them
 * in place, which reads xyz twice and needs no extra memory. After
 * seams_system_stage_positions(sys, 1) it reads xyz once instead (see
 * there).
 *
 * @param sys Frame handle with a stored frame.
 * @param xyz n times 3 Cartesian coordinates, Angstrom. Every value must
 *        be finite.
 * @param n Atom count. Must equal the stored count.
 * @param box9 New cell (ox, oy, oz, lx, ly, lz, xy, xz, yz), or NULL to
 *        keep the stored cell.
 * @param box_len 9 when box9 is not NULL; ignored otherwise.
 * @return ::SEAMS_OK; ::SEAMS_ERR_NULL for a NULL sys or xyz;
 *         ::SEAMS_ERR_STATE when no frame is stored; ::SEAMS_ERR_RANGE
 *         when n differs from the stored count or box_len != 9;
 *         ::SEAMS_ERR_BOX for a non-finite cell or a non-positive edge;
 *         ::SEAMS_ERR_ARG for a non-finite coordinate; ::SEAMS_ERR_NOMEM;
 *         ::SEAMS_ERR_INTERNAL.
 * @par Handle access
 * Non-const: replaces the stored coordinates and cell.
 * @par Ownership
 * The library copies xyz and box9. The caller keeps them and may reuse
 * the buffers as soon as the call returns.
 * @par Thread safety
 * Exclusive use of sys.
 * @par Example
 * @code
 * seams_system_set_atoms(sys, frame0, n, box9, 9, types, NULL);
 * for (size_t f = 1; f < nframes; ++f) {
 *   seams_system_set_positions(sys, frames[f], n, NULL, 0);
 *   seams_chill_plus_n(sys, 3.5, 0, SEAMS_GRAPH_CUTOFF, SEAMS_CHILL_PUBLISHED,
 *                      labels, n, NULL);
 * }
 * @endcode
 */
SEAMS_API seams_status seams_system_set_positions(seams_system *sys, const double *xyz,
                                                  size_t n, const double *box9,
                                                  size_t box_len);

/**
 * @brief Trade memory for a single pass in seams_system_set_positions().
 *
 * With enable nonzero the handle keeps a second set of per-atom records,
 * about 72 bytes per atom. seams_system_set_positions() then checks and
 * writes the coordinates into that set in one pass and swaps it with the
 * stored frame when every value was finite, so a failed call still leaves
 * the frame unchanged. The first reload after each load fills the second
 * set from the frame. With enable 0, the default, the second set is freed.
 * The setting survives loads.
 *
 * @param sys Frame handle, with or without a stored frame.
 * @param enable Nonzero to keep the second set, 0 to free it.
 * @return ::SEAMS_OK; ::SEAMS_ERR_NULL for a NULL sys.
 * @par Handle access
 * Non-const: changes how later reloads use memory. The stored frame is
 * not touched.
 * @par Ownership
 * The library owns the second set and frees it with enable 0 or
 * seams_system_destroy().
 * @par Thread safety
 * Exclusive use of sys.
 * @par Example
 * @code
 * seams_system_stage_positions(sys, 1);
 * for (size_t f = 1; f < nframes; ++f) {
 *   seams_system_set_positions(sys, frames[f], n, NULL, 0);
 * }
 * @endcode
 */
SEAMS_API seams_status seams_system_stage_positions(seams_system *sys, int enable);

/**
 * @brief Copy the stored cell out as nine LAMMPS numbers.
 * @param sys Frame handle.
 * @param box9 Caller buffer: ox, oy, oz, lx, ly, lz, xy, xz, yz. The
 *        values match the cell that was stored, to rounding.
 * @param box_len Must be 9.
 * @return ::SEAMS_OK; ::SEAMS_ERR_NULL; ::SEAMS_ERR_RANGE when
 *         box_len != 9; ::SEAMS_ERR_STATE when no frame is stored.
 * @par Ownership
 * The caller owns box9.
 * @par Thread safety
 * Read-only: may run alongside other const calls on sys, and not
 * alongside a non-const call on sys.
 * @par Example
 * @code
 * double cell[9];
 * seams_system_box(sys, cell, 9);
 * @endcode
 */
SEAMS_API seams_status seams_system_box(const seams_system *sys, double *box9,
                                        size_t box_len);

/**
 * @brief Copy the stored coordinates out.
 * @param sys Frame handle.
 * @param xyz Caller buffer for 3 n doubles, Angstrom, in stored order.
 * @param len Length of xyz in doubles. Must be at least 3 n.
 * @return ::SEAMS_OK; ::SEAMS_ERR_NULL; ::SEAMS_ERR_STATE when no frame
 *         is stored; ::SEAMS_ERR_RANGE for a short buffer.
 * @par Ownership
 * The caller owns xyz.
 * @par Thread safety
 * Read-only: may run alongside other const calls on sys, and not
 * alongside a non-const call on sys.
 * @par Example
 * @code
 * size_t n = 0;
 * seams_system_natoms(sys, &n);
 * double *xyz = malloc(3 * n * sizeof *xyz);
 * seams_system_positions(sys, xyz, 3 * n);
 * @endcode
 */
SEAMS_API seams_status seams_system_positions(const seams_system *sys, double *xyz,
                                              size_t len);

/**
 * @brief Copy the stored atom types out.
 * @param sys Frame handle.
 * @param types Caller buffer for n ints, in stored order.
 * @param len Length of types. Must be at least n.
 * @return ::SEAMS_OK; ::SEAMS_ERR_NULL; ::SEAMS_ERR_STATE when no frame
 *         is stored; ::SEAMS_ERR_RANGE for a short buffer.
 * @par Ownership
 * The caller owns types.
 * @par Thread safety
 * Read-only: may run alongside other const calls on sys, and not
 * alongside a non-const call on sys.
 * @par Example
 * @code
 * seams_system_types(sys, types, n);
 * @endcode
 */
SEAMS_API seams_status seams_system_types(const seams_system *sys, int *types, size_t len);

/**
 * @brief Copy the stored atom ids out.
 * @param sys Frame handle.
 * @param ids Caller buffer for n ints, in stored order.
 * @param len Length of ids. Must be at least n.
 * @return ::SEAMS_OK; ::SEAMS_ERR_NULL; ::SEAMS_ERR_STATE when no frame
 *         is stored; ::SEAMS_ERR_RANGE for a short buffer.
 * @par Ownership
 * The caller owns ids.
 * @par Thread safety
 * Read-only: may run alongside other const calls on sys, and not
 * alongside a non-const call on sys.
 * @par Example
 * @code
 * seams_system_ids(sys, ids, n);
 * @endcode
 */
SEAMS_API seams_status seams_system_ids(const seams_system *sys, int *ids, size_t len);

/** @} */

#ifdef __cplusplus
}
#endif

#endif
