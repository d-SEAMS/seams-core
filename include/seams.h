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
#include <stdint.h>

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
  SEAMS_ERR_RANGE = 2,     /**< A number or a buffer length is outside the accepted range, a tilted cell needs a lattice coefficient at or above 2^52 for its minimum image, or the cell is outside what the linked-cell search resolves, such as an edge that needs more bins than linkcell allows. */
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

/** @brief Neighbour convention for tetrahedral order. */
typedef enum seams_tetrahedral_neighbours {
  SEAMS_TETRAHEDRAL_NEAREST_FOUR = 0, /**< Four nearest selected partners; cutoff is ignored. */
  SEAMS_TETRAHEDRAL_WITHIN_CUTOFF = 1 /**< Every selected partner at or inside cutoff. */
} seams_tetrahedral_neighbours;

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

/**
 * @brief Collective variable evaluated by seams_cv_value().
 *
 * Ring sizes in arg are 0 (every size) or 3 through the workspace's
 * max_ring. The ring kinds and the hexagon pairs are polynomials in the
 * bond weights, as smooth as the switch of seams_switch_rational(). The
 * gyration kinds are smooth only for a compact selection of positive
 * total weight, and the volume is in general not differentiable where the
 * cluster becomes flat; seams_cv_value() gives the conditions. A state is
 * as smooth as the variables its clauses read.
 */
typedef enum seams_cv_kind {
  SEAMS_CV_RING_COUNT = 0,      /**< Sum of soft ring weights of size arg; as smooth as the switch. */
  SEAMS_CV_RING_ATOMS = 1,      /**< Soft number of selected atoms on a ring of size arg; as smooth as the switch. */
  SEAMS_CV_HEX_BASAL_PAIRS = 2, /**< Soft number of hexagon pairs that pass the HC basal test, arg 0; as smooth as the switch. */
  SEAMS_CV_GYRATION_RG2 = 3,    /**< Weighted squared radius of gyration, square angstrom; smooth for a compact selection of positive weight. */
  SEAMS_CV_GYRATION_VOLUME = 4, /**< Volume of the gyration ellipsoid, cubic angstrom; as RG2, and in general not differentiable where the cluster becomes flat. */
  SEAMS_CV_STATE = 5            /**< Soft indicator of state arg from seams_cv_state_clause(); as smooth as its clauses. */
} seams_cv_kind;

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
/** @brief Streaming mean-squared displacement accumulator. Free with seams_msd_destroy(). */
typedef struct seams_msd seams_msd;
/** @brief Streaming hydrogen-bond lifetime accumulator. Free with seams_hbond_tcf_destroy(). */
typedef struct seams_hbond_tcf seams_hbond_tcf;
/** @brief Streaming orientational time-correlation accumulator. Free with seams_reorient_destroy(). */
typedef struct seams_reorient seams_reorient;
/** @brief Smooth collective-variable workspace. Free with seams_cv_destroy(). */
typedef struct seams_cv seams_cv;

/**
 * @brief Callback that returns quenched coordinates for one frame.
 *
 * The callback receives one frame and writes the quenched coordinates into
 * @p xyz_out (3 * n doubles) and, optionally, a replacement box into
 * @p box9_out (nine doubles). Returning ::SEAMS_OK accepts @p xyz_out as the
 * quenched frame. Returning an error aborts the quenched analysis call.
 *
 * @param xyz Input coordinates, 3 * n doubles.
 * @param n Atom count.
 * @param box9 Input cell, nine doubles (ox, oy, oz, lx, ly, lz, xy, xz, yz).
 * @param box_len Always 9.
 * @param xyz_out Output coordinates, 3 * n doubles.
 * @param xyz_out_len Always 3 * n.
 * @param box9_out Optional output box. May be NULL to keep @p box9.
 * @param box9_out_len 9 when @p box9_out is non-NULL, else 0.
 * @param user_data Opaque caller pointer passed through unchanged.
 * @return ::SEAMS_OK, or another ::seams_status to reject this frame.
 */
typedef seams_status (*seams_quench_callback)(const double *xyz, size_t n, const double *box9,
                                              size_t box_len, double *xyz_out, size_t xyz_out_len,
                                              double *box9_out, size_t box9_out_len,
                                              void *user_data);

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
 * @brief Version text of the built library.
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
 *
 * Type 0 on the atoms (1, 1, 1) and (2, 1.2, 1.1) in an orthorhombic
 * cell of edge 1e160 with origin (0, 0, 0) returns ::SEAMS_OK from the
 * cutoff graph: an empty list at cutoff 0.5 and each other at cutoff 3.5,
 * the same rows as that pair and type on an edge of 30 with that origin.
 * The same atoms and type with origin (5, -2, 0.25) list each other at
 * cutoff 0.5 on the edge of 1e160 and do not on an edge of 30. Type 0 on
 * (-0.9, 0, 0) with (0.1, 0.2, 0.1), and on (-2, -1, -4) with
 * (-1, -0.8, -3.9), returns ::SEAMS_OK and both neighbours at cutoff 0.5
 * on edge 1e160 with origin (0, 0, 0), and an empty list on edge 30.
 * Type 0 on (0, 0, 0) and 0x1.6c2d4256ffcc4p+588 returns ::SEAMS_OK and
 * both neighbours at cutoff 0.5 on that edge and origin.
 *
 * When the orthorhombic cutoff graph bins a coordinate, the fraction is
 * (coordinate - origin) times the binary64 reciprocal of the edge. When
 * that fraction is outside [0, 1) it subtracts the floor. A result that
 * is not strictly below 1, NaN included, becomes 0. The stored coordinate
 * is the origin plus the fraction times the edge. Type 0 on (0, 0, 0) and
 * (1e162, 0, 0), origin (0, 0, 0), edge 1e160, cutoff 0.5, returns
 * ::SEAMS_OK and an empty list. The product is 99.99999999999999. The
 * floor leaves 0.9999999999999858, and the stored x coordinate is
 * 9.999999999999858e159. Type 0 on (0, 0, 0) and (1.5e160, 0, 0), that
 * edge and origin, cutoff 0.5, returns ::SEAMS_OK and an empty list. The
 * product is 1.5. The floor leaves 0.5, and the stored x coordinate is
 * 5e159. Type 0 on (3.2581295480326344e+46, 0, 0) and
 * (3.2581295480326338e+46, 0, 0), that edge and origin, cutoff 0.5,
 * returns ::SEAMS_OK and both neighbours. Both stored coordinates are
 * 3.2581295480326338e+46. Type 0 on (5, 0, 0) and (4096, 0, 0), origin
 * (5, 0, 0), edge 1e160, cutoff 4091, returns ::SEAMS_OK and both
 * neighbours. The wrap stores x = 4096 at 4095.9999999999995. Type 0 on
 * (-1e308, 0, 0) and (1e308, 0, 0), origin (-1e308, 0, 0), edge 30,
 * cutoff 0.5, returns ::SEAMS_OK and both neighbours. The product for
 * the second atom is inf, inf minus its floor is NaN, and that NaN
 * becomes 0, so both atoms are stored at the origin. Type 0 on (0, 0, 0)
 * and (1, 0, 0), origin (0, 0, 0), edge 2^-1070, cutoff 0.5, returns
 * ::SEAMS_OK and both neighbours. The reciprocal of that edge is inf,
 * both fractions become 0, and both atoms are stored at the origin.
 *
 * Type 0 on the mutual graph returns ::SEAMS_ERR_RANGE for the atoms
 * (1, 1, 1) and (2, 1.2, 1.1) on the orthorhombic cell of edge 1e160 with
 * origin (0, 0, 0), and on the triclinic cell of that edge with that
 * origin and tilts 1, 0.5 and -0.25, at cutoffs 0.5 and 3.5. Type 99 on
 * each of those frames returns ::SEAMS_OK and an empty list at those
 * cutoffs. On that triclinic cell, type 0 on the cutoff graph returns
 * ::SEAMS_OK and both neighbours for those atoms at cutoff 3.5, and
 * ::SEAMS_ERR_RANGE at cutoff 0.5. Type 99 at that cutoff returns
 * ::SEAMS_OK and an empty list.
 *
 * Type 0 on (0, 0, 0) and (10 * 2^52, 0, 0), edge 10, origin (0, 0, 0)
 * and tilt xy = 1, returns ::SEAMS_ERR_RANGE from the cutoff graph at
 * cutoffs 0.5 and 3.5. Type 0 on (5, 5, 5) and (45035996273704968, 5, 5)
 * on that cell returns ::SEAMS_OK and an empty list at cutoff 0.5, and
 * ::SEAMS_ERR_RANGE at cutoff 3.5. Type 0 on (0, 0, 0) and (2^53 - 1, 0, 0),
 * edges 2, 2 and 2, origin (0, 0, 0) and tilt xy = 0.25, returns
 * ::SEAMS_ERR_RANGE at cutoffs 0.5 and 1.0, with the message that the
 * minimum-image wrap is not a lattice image, and ::SEAMS_OK with both
 * neighbours at cutoff 3.5. Type 0 on (0, 0, 0) and
 * (50 * 2^52, 10 * 2^51, 0), edges 10, 10 and 10, origin (0, 0, 0) and
 * tilt xy = 100, returns ::SEAMS_OK and both neighbours at cutoffs 0.5
 * and 3.5. Type 0 on (0, 0, 0) and (10 * 2^52, 0, 0) in an orthorhombic
 * cell of edge 10 with origin (0, 0, 0) returns ::SEAMS_OK and both
 * neighbours at cutoffs 0.5 and 3.5. Type 0 on (4, 0, 0) and
 * (4 + 10 * 2^52, 0, 0) in that cell returns ::SEAMS_OK and an empty list
 * at both cutoffs. Type 0 on (0, 0, 0) and (0.5, 0, 0) in that cell
 * returns ::SEAMS_OK and an empty list at cutoff 0.5.
 *
 * For a pair the orthorhombic cutoff cell list visits, it keeps the
 * neighbour only when the squared separation of the wrapped coordinates
 * and a periodic cell shift is strictly below the squared cutoff, both
 * computed in binary64. An overflowing square and a NaN square fail that
 * test. Type 0 on (0, 0, 0) and (2e154, 0, 0), edge 1e160, origin
 * (0, 0, 0), cutoff 3e154, returns ::SEAMS_OK and an empty list. Type 0
 * on (0.9e308, 0, 0) and (0.901e308, 0, 0), origin (1e308, 0, 0), edge
 * 1e308, cutoff 0.5, returns ::SEAMS_OK and an empty list. Both stored x
 * coordinates are inf, so each squared separation is NaN.
 *
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
 * free(off);
 * free(cols);
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

/**
 * @brief Errington-Debenedetti tetrahedral order for arbitrary selections.
 *
 * centres and partners contain zero-based frame indices. A NULL pointer with
 * zero length selects every atom. A centre is excluded from its own partner
 * shell. q = 1 - 9/(2 N (N-1)) sum_{j<k} (cos psi_jk + 1/3)^2 over the N
 * shell neighbours: a perfect tetrahedron is 1 and an ideal gas averages 0.
 * NEAREST_FOUR uses the four nearest selected partners, N = 4 and the
 * Errington-Debenedetti factor 3/8, and ignores cutoff. WITHIN_CUTOFF uses
 * every selected partner at or inside cutoff. In both modes a centre with
 * fewer than four neighbours receives quiet NaN.
 *
 * @param sys Frame handle.
 * @param centres Centre indices, or NULL when n_centres is zero.
 * @param n_centres Centre count. Zero means every atom.
 * @param partners Partner indices, or NULL when n_partners is zero.
 * @param n_partners Partner count. Zero means every atom.
 * @param convention One ::seams_tetrahedral_neighbours value. Any other
 *        value returns ::SEAMS_ERR_ARG.
 * @param cutoff Positive finite cutoff in WITHIN_CUTOFF mode; ignored in
 *        NEAREST_FOUR mode.
 * @param q Caller buffer, one value per selected centre.
 * @param q_len Length of q.
 * @return ::SEAMS_OK, or ::SEAMS_ERR_NULL, ::SEAMS_ERR_STATE,
 *         ::SEAMS_ERR_RANGE, ::SEAMS_ERR_ARG, ::SEAMS_ERR_NOMEM,
 *         ::SEAMS_ERR_INTERNAL.
 * @par Handle access
 * Non-const only because the call clears the handle error slot.
 * @par Ownership
 * The caller owns every selection and output buffer.
 * @par Thread safety
 * Exclusive use of sys.
 */
SEAMS_API seams_status seams_tetrahedral_order(
    seams_system *sys, const int *centres, size_t n_centres,
    const int *partners, size_t n_partners, int convention, double cutoff,
    double *q, size_t q_len);

/**
 * @brief Local Structure Index for arbitrary selections.
 *
 * Shiratani-Sasai index. For each centre, selected-partner distances are
 * sorted, r_1 < ... < r_N <= cutoff < r_{N+1}. The result is the population
 * variance of the N gaps r_{i+1} - r_i, i = 1..N; the last gap reaches the
 * first partner beyond cutoff. A centre with no partner inside cutoff, or
 * none beyond it, receives quiet NaN.
 *
 * @param sys Frame handle.
 * @param centres Centre indices, or NULL when n_centres is zero.
 * @param n_centres Centre count. Zero means every atom.
 * @param partners Partner indices, or NULL when n_partners is zero.
 * @param n_partners Partner count. Zero means every atom.
 * @param cutoff Positive finite shell cutoff.
 * @param lsi Caller buffer, one value per selected centre.
 * @param lsi_len Length of lsi.
 * @return ::SEAMS_OK, or ::SEAMS_ERR_NULL, ::SEAMS_ERR_STATE,
 *         ::SEAMS_ERR_RANGE, ::SEAMS_ERR_ARG, ::SEAMS_ERR_NOMEM,
 *         ::SEAMS_ERR_INTERNAL.
 * @par Handle access
 * Non-const only because the call clears the handle error slot.
 * @par Ownership
 * The caller owns every selection and output buffer.
 * @par Thread safety
 * Exclusive use of sys.
 */
SEAMS_API seams_status seams_lsi(seams_system *sys, const int *centres,
                                 size_t n_centres, const int *partners,
                                 size_t n_partners, double cutoff, double *lsi,
                                 size_t lsi_len);

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
 * @brief Primitive rings from an explicit caller-supplied CSR graph.
 *
 * offsets has n_vertices + 1 entries and starts at zero. columns contains
 * zero-based vertex indices. Directed input is symmetrized by edge union;
 * self edges and duplicates are ignored. The Franzblau search is the same
 * implementation used by seams_rings_find(). geometry may be NULL. When
 * supplied, its atom count must equal n_vertices and winding queries are
 * enabled.
 *
 * @param offsets CSR offsets.
 * @param offsets_len Number of offsets; at least one.
 * @param columns CSR columns, or NULL when columns_len is zero.
 * @param columns_len Number of columns.
 * @param max_depth Largest ring size, from 3 through 12.
 * @param geometry Optional frame supplying coordinates and periodic cell.
 * @param out Receives a new ring handle.
 * @return ::SEAMS_OK, or ::SEAMS_ERR_NULL, ::SEAMS_ERR_STATE,
 *         ::SEAMS_ERR_RANGE, ::SEAMS_ERR_ARG, ::SEAMS_ERR_NOMEM,
 *         ::SEAMS_ERR_INTERNAL.
 * @par Handle access
 * Read-only access to geometry.
 * @par Ownership
 * The caller owns all arrays and frees *out with seams_rings_destroy().
 * @par Thread safety
 * Read-only calls may share geometry. The returned handle is independent.
 */
SEAMS_API seams_status seams_rings_from_csr(
    const int *offsets, size_t offsets_len, const int *columns,
    size_t columns_len, int max_depth, const seams_system *geometry,
    seams_rings **out);

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
 * @return ::SEAMS_OK, or ::SEAMS_ERR_NULL, ::SEAMS_ERR_RANGE,
 *         ::SEAMS_ERR_NOMEM, ::SEAMS_ERR_INTERNAL.
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
 * @brief Periodic-axis winding and contractibility of every ring.
 *
 * axis_winding receives the oriented integer winding along axis (0, 1, or
 * 2): the minimum-image bond steps summed in the order seams_rings_vertices()
 * lists the ring, in units of that cell vector. A ring that crosses the cell
 * once along +x is +1 on axis 0; the reverse traversal is -1. contractible
 * receives 1 exactly when all three winding components are zero. Either
 * output may be NULL, but not both. Geometry must
 * have been supplied to seams_rings_from_csr(), or carried by
 * seams_rings_find().
 * @param rings Ring handle.
 * @param axis Periodic cell axis, 0, 1, or 2.
 * @param axis_winding Signed winding output, or NULL.
 * @param axis_winding_len Length of axis_winding.
 * @param contractible One for a zero three-component winding, or NULL.
 * @param contractible_len Length of contractible.
 * @return ::SEAMS_OK, or ::SEAMS_ERR_NULL, ::SEAMS_ERR_STATE,
 *         ::SEAMS_ERR_RANGE.
 * @par Ownership
 * The caller owns both output buffers.
 * @par Thread safety
 * Read-only: may run alongside other const calls on rings.
 */
SEAMS_API seams_status seams_rings_winding(
    const seams_rings *rings, int axis, int *axis_winding,
    size_t axis_winding_len, int *contractible, size_t contractible_len);

/**
 * @brief Chordless cycles of one size.
 *
 * A cycle is chordless when no edge of the graph joins two of its
 * vertices that are not next to each other on the cycle. Every primitive
 * ring is chordless; a chordless cycle need not be primitive. Each cycle
 * is written once as size atom indices: it starts at its lowest index
 * and runs towards the smaller of that atom's two cycle neighbours.
 * Cycles are sorted lexicographically, so the list does not depend on
 * the thread count.
 * @param nb Graph from seams_neighbours_build().
 * @param size Cycle length, from 3 through 12.
 * @param rings Caller buffer of size times the cycle count ints, or
 *        NULL to skip the list.
 * @param rings_len Length of rings in ints. Ignored when rings is NULL.
 * @param n_rings Receives the cycle count. Required when rings is not
 *        NULL.
 * @param per_atom Caller buffer of one int per atom, or NULL. Receives
 *        the number of cycles through each atom.
 * @param per_atom_len Length of per_atom. Ignored when per_atom is NULL.
 * @return ::SEAMS_OK, or ::SEAMS_ERR_NULL, ::SEAMS_ERR_RANGE,
 *         ::SEAMS_ERR_NOMEM, ::SEAMS_ERR_INTERNAL. rings, n_rings and
 *         per_atom may not all be NULL. A short rings buffer returns
 *         ::SEAMS_ERR_RANGE after n_rings and per_atom are written.
 * @par Ownership
 * The caller owns rings and per_atom. Indices are positions in the frame.
 * @par Thread safety
 * Exclusive use of nb.
 * @par Example
 * @code
 * size_t n_six = 0;
 * seams_rings_chordless(nb, 6, NULL, 0, &n_six, NULL, 0);
 * @endcode
 */
SEAMS_API seams_status seams_rings_chordless(const seams_neighbours *nb, int size, int *rings,
                                             size_t rings_len, size_t *n_rings,
                                             int *per_atom, size_t per_atom_len);

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
 * @defgroup capi_cv Smooth collective variables
 * @brief Soft ring counts, hexagon-pair counts, gyration and states, with gradients.
 *
 * Every bond within d_max carries the weight s(r) of seams_switch_rational().
 * A candidate ring is any simple cycle of at most max_ring members in the
 * graph at d_max + skin. Its weight is the product of its bond weights times
 * (1 - W_P) for every shortcut path P between two of its members, W_P the
 * product of the path's bond weights. With bond weights 0 and 1 the rings of
 * weight one are the primitive rings of seams_rings_find(). The ring counts,
 * ring-atom counts and hexagon pairs are fixed polynomials in the pair
 * weights, so they are as smooth as s(r) also where a bond crosses the list
 * cutoff or a ring gains or loses a shortcut.
 * The weights are twice continuously differentiable in the coordinates for
 * n >= 3, and for n = 2 with d0 = 0. At d0 > 0 the weight is only
 * C^(n-1) at r = d0; with d0 = 0 and n = 1 it has a cone, s = 1 - r / r0
 * to first order, where two atoms coincide.
 *
 * The gyration values are not polynomials in the weights. They are as
 * smooth as the weights for a compact selection, one whose atoms are all
 * closer to each other than half the narrowest cell width in some periodic
 * image, while the total weight is positive. Otherwise the unwrapped
 * positions can jump by a cell vector, and with ring weights the value
 * drops to 0 where the total weight reaches 0. The volume is a square root
 * of det G: continuous where G is, but in general not differentiable where
 * the cluster becomes flat, det G = 0. It has a kink there, like the
 * distance of one atom from the plane of the others.
 *
 * Gradients are d value / d x, 3 * natoms doubles in frame order,
 * x0 y0 z0 x1 ..., zero for atoms outside the selection. The force of a
 * bias V(value) on atom i is -dV/dvalue times its three entries.
 * @{
 */

/**
 * @brief PLUMED RATIONAL switching function with a C2 taper to zero.
 *
 * For x = (r - d0) / r0, s = (1 - x^n) / (1 - x^m) for r > d0 and 1 for
 * r <= d0, times a taper that is 1 below d_on, 0 at and beyond d_max, and
 * 1 - u^3 (10 - 15u + 6u^2) for u = (r - d_on) / (d_max - d_on) between
 * them. s and its first two derivatives vanish at d_max. On r >= 0, s is
 * C^2 in r except at r = d0 > 0 with n < 3, where it is only C^(n-1). At
 * r = d0, dsdr is the right derivative: -1 / r0 for n = 1 and 0 otherwise.
 * Above x = 1, s is evaluated in 1 / x, so for every finite argument s is
 * finite when n and m are at most 4096. A larger exponent is
 * ::SEAMS_ERR_RANGE and the message names 4096. ds/dr can exceed double range when r0 or
 * d_max - d_on is near the smallest doubles; when dsdr is not NULL that
 * is ::SEAMS_ERR_RANGE, and neither output is written.
 * @param r Distance, Angstrom. Not negative.
 * @param r0 Positive length scale, Angstrom.
 * @param d0 Plateau radius, Angstrom, at least zero.
 * @param n Numerator exponent, from 1 through 4096.
 * @param m Denominator exponent, greater than n and at most 4096. Zero
 *        selects 2 * n, which needs 2 * n <= 4096.
 * @param d_on Taper start, Angstrom, from d0 and below d_max.
 * @param d_max Radius where s reaches zero, Angstrom.
 * @param s Receives s(r).
 * @param dsdr Receives ds/dr in 1/Angstrom, the right derivative at
 *        r = d0, or NULL.
 * @return ::SEAMS_OK, or ::SEAMS_ERR_NULL, ::SEAMS_ERR_RANGE.
 * @par Ownership
 * The caller owns s and dsdr. The call keeps no state.
 * @par Thread safety
 * Safe from any thread.
 * @par Example
 * @code
 * double sw = 0, dsw = 0;
 * seams_switch_rational(3.0, 0.5, 2.5, 6, 0, 3.3, 3.5, &sw, &dsw);
 * @endcode
 */
SEAMS_API seams_status seams_switch_rational(double r, double r0, double d0, int n, int m,
                                             double d_on, double d_max, double *s,
                                             double *dsdr);

/**
 * @brief Create a smooth collective-variable workspace.
 * @param r0 Switching length scale, Angstrom. See seams_switch_rational().
 * @param d0 Switching plateau radius, Angstrom.
 * @param n Switching numerator exponent, from 1 through 4096.
 * @param m Switching denominator exponent, greater than n and at most 4096.
 *        Zero selects 2 * n, which needs 2 * n <= 4096.
 * @param d_on Taper start, Angstrom.
 * @param d_max Radius where every bond weight reaches zero, Angstrom.
 * @param max_ring Largest candidate ring, from 3 through 10.
 * @param out Receives a new handle. On failure *out is NULL.
 * @return ::SEAMS_OK, or ::SEAMS_ERR_NULL, ::SEAMS_ERR_RANGE,
 *         ::SEAMS_ERR_NOMEM.
 * @par Ownership
 * The caller frees *out with seams_cv_destroy(). The handle holds no
 * frame until seams_cv_rebuild().
 * @par Thread safety
 * The new handle is private to the caller.
 * @par Example
 * @code
 * seams_cv *cv = NULL;
 * if (seams_cv_create(0.5, 2.5, 6, 0, 3.3, 3.5, 6, &cv) == SEAMS_OK) {
 *   seams_cv_destroy(cv);
 * }
 * @endcode
 */
SEAMS_API seams_status seams_cv_create(double r0, double d0, int n, int m, double d_on,
                                       double d_max, int max_ring, seams_cv **out);

/**
 * @brief Free a collective-variable workspace.
 * @param cv Handle from seams_cv_create(), or NULL. NULL is ignored.
 * @par Ownership
 * Releases the handle, its topology and its state clauses.
 * @par Thread safety
 * The caller must not use cv after this returns.
 * @par Example
 * @code
 * seams_cv *cv = NULL;
 * seams_cv_create(0.3, 2.8, 6, 0, 2.9, 3.0, 6, &cv);
 * seams_cv_destroy(cv);
 * @endcode
 */
SEAMS_API void seams_cv_destroy(seams_cv *cv);

/**
 * @brief List the graph and fix the candidate rings for the next frames.
 *
 * Bonds are the selected pairs within d_max + skin, from linkcell. The
 * candidate rings, their shortcut paths and the hexagon pairs are fixed
 * until the next rebuild. A ring, ring-atom or hexagon-pair value computed
 * from this topology equals the value from a fresh rebuild while no
 * selected atom has moved more than skin / 2 from its position here, the
 * cell is unchanged, and every atom type is unchanged; so does a gyration
 * value of a compact selection (see seams_cv_value()). A type change marks
 * the selection stale. The ring and shortcut searches are bounded by the
 * ring size. A graph that would pass 20000000 ring-search steps, 2000000
 * candidate rings, 20000000 shortcut-search steps, or 2000000 shortcut
 * paths is ::SEAMS_ERR_RANGE and the message names that limit; rings and
 * paths are not dropped. seams_cv_needs_rebuild() reports when the
 * topology no longer matches the frame, and seams_cv_value() then refuses
 * the frame.
 * @param cv Workspace handle.
 * @param sys Frame handle that holds a frame.
 * @param type Species to select. Zero and negative select every atom.
 * @param atoms Frame indices to select, or NULL for every atom of type.
 * @param n_atoms Length of atoms. Ignored when atoms is NULL.
 * @param skin Verlet skin, Angstrom, at least zero. d_max + skin must be
 *        below half the narrowest width of the cell.
 * @return ::SEAMS_OK, or ::SEAMS_ERR_NULL, ::SEAMS_ERR_STATE,
 *         ::SEAMS_ERR_RANGE, ::SEAMS_ERR_ARG, ::SEAMS_ERR_BOX,
 *         ::SEAMS_ERR_NOMEM, ::SEAMS_ERR_INTERNAL.
 * @par Handle access
 * Writes the topology into cv. Reads sys.
 * @par Ownership
 * The caller owns atoms. The workspace keeps a copy of the selection and of
 * the coordinates and cell of this frame.
 * @par Thread safety
 * Exclusive use of cv. Read-only use of sys.
 * @par Example
 * @code
 * seams_cv *cv = NULL;
 * seams_cv_create(0.3, 2.8, 6, 0, 2.9, 3.0, 6, &cv);
 * seams_cv_rebuild(cv, sys, 0, NULL, 0, 0.1);
 * seams_cv_destroy(cv);
 * @endcode
 */
SEAMS_API seams_status seams_cv_rebuild(seams_cv *cv, const seams_system *sys, int type,
                                        const int *atoms, size_t n_atoms, double skin);

/**
 * @brief Whether a frame needs seams_cv_rebuild() before the next value.
 *
 * It compares state only and does not bin the cell. A too-fine mesh is
 * reported by seams_cv_rebuild() as ::SEAMS_ERR_RANGE. stale is written
 * only on ::SEAMS_OK and left unchanged on any error. With skin 0.1, the
 * atoms (0, 0, 0) and (1, 0, 0) in an orthorhombic cell of edge 30 with
 * origin (0, 0, 0) rebuild, and moving the first atom to (30 * 2^52, 0, 0)
 * on that same cell returns ::SEAMS_ERR_RANGE and leaves a pre-set stale
 * of 7 unchanged.
 * @param cv Workspace handle.
 * @param sys Frame handle that holds a frame.
 * @param stale Written only on ::SEAMS_OK. Then 1 when the workspace has
 *        not been rebuilt, the atom count, cell, or origin differs, an
 *        atom type differs, or a selected atom's squared minimum-image
 *        displacement is greater than (skin / 2)^2, and 0 when none of
 *        those holds. Left unchanged on any error.
 * @return ::SEAMS_OK, or ::SEAMS_ERR_NULL, ::SEAMS_ERR_STATE,
 *         ::SEAMS_ERR_RANGE.
 * @par Ownership
 * The caller owns stale.
 * @par Thread safety
 * Read-only: may run alongside other const calls on cv and sys.
 * @par Example
 * @code
 * seams_cv *cv = NULL;
 * seams_cv_create(0.3, 2.8, 6, 0, 2.9, 3.0, 6, &cv);
 * seams_cv_rebuild(cv, sys, 0, NULL, 0, 0.1);
 * int stale = 0;
 * seams_cv_needs_rebuild(cv, sys, &stale);
 * seams_cv_destroy(cv);
 * @endcode
 */
SEAMS_API seams_status seams_cv_needs_rebuild(const seams_cv *cv, const seams_system *sys,
                                              int *stale);

/**
 * @brief Sizes of the fixed topology.
 * @param cv Workspace handle.
 * @param n_edges Receives the number of listed bonds, or NULL.
 * @param n_rings Receives the number of candidate rings, or NULL.
 * @param n_pairs Receives the number of candidate hexagon pairs, or NULL.
 * @return ::SEAMS_OK, or ::SEAMS_ERR_NULL, ::SEAMS_ERR_STATE before the
 *         first rebuild.
 * @par Ownership
 * The caller owns the three outputs.
 * @par Thread safety
 * Read-only: may run alongside other const calls on cv.
 * @par Example
 * @code
 * seams_cv *cv = NULL;
 * seams_cv_create(0.3, 2.8, 6, 0, 2.9, 3.0, 6, &cv);
 * seams_cv_rebuild(cv, sys, 0, NULL, 0, 0.1);
 * size_t edges = 0, rings = 0, pairs = 0;
 * seams_cv_topology(cv, &edges, &rings, &pairs);
 * seams_cv_destroy(cv);
 * @endcode
 */
SEAMS_API seams_status seams_cv_topology(const seams_cv *cv, size_t *n_edges,
                                         size_t *n_rings, size_t *n_pairs);

/**
 * @brief Value and gradient of one collective variable.
 *
 * ::SEAMS_CV_RING_COUNT sums the ring weights of size arg.
 * ::SEAMS_CV_RING_ATOMS sums 1 - prod (1 - w) over the rings of size arg
 * through each selected atom. ::SEAMS_CV_HEX_BASAL_PAIRS sums
 * w_a w_b B over vertex-disjoint hexagon pairs, B the soft form of the
 * hexagonal-cage basal test of either ring against the other; arg is 0.
 * These three are as smooth as the switch.
 * ::SEAMS_CV_GYRATION_RG2 and ::SEAMS_CV_GYRATION_VOLUME use the gyration
 * tensor G of the selected atoms, unwrapped along bonds, with weight 1
 * when arg is -1 and the ring participation of size arg otherwise; the
 * volume is (4 pi / 3) 5^(3/2) sqrt(det G), the volume of the uniform
 * ellipsoid with that tensor. A part of the selection that no bond reaches
 * is placed by its minimum image from the first selected atom. Both are
 * smooth while the selection is compact, its atoms all closer to each
 * other than half the narrowest cell width in some periodic image, and its
 * total weight is positive. With total weight 0 both values and gradients
 * are 0, so with ring weights they jump to 0 as the last ring breaks. The
 * volume is in general not differentiable where the cluster becomes flat,
 * det G = 0; there it and its gradient are 0. The same 0 value and 0
 * gradient are returned when the smallest eigenvalue of G is at most
 * 1e-10 times the largest, so a coplanar cloud whose determinant rounds
 * positive stays flat. ::SEAMS_CV_STATE returns
 * state arg, as smooth as the variables its clauses read.
 * A value or gradient entry outside double range, from switch parameters
 * near the smallest doubles or coordinates near the largest, is
 * ::SEAMS_ERR_RANGE, and neither value nor grad is written.
 * @param cv Workspace handle after seams_cv_rebuild().
 * @param sys Frame with the atom count and cell of the rebuild.
 * @param kind A ::seams_cv_kind code.
 * @param arg Ring size, weight choice, or state index, as above.
 * @param value Receives the value.
 * @param grad Caller buffer of 3 * natoms doubles for d value / d x in
 *        1/Angstrom times the unit of value, or NULL to skip it.
 * @param grad_len Length of grad. Ignored when grad is NULL.
 * @return ::SEAMS_OK, or ::SEAMS_ERR_NULL, ::SEAMS_ERR_STATE when cv
 *         needs a rebuild for sys, ::SEAMS_ERR_RANGE, ::SEAMS_ERR_NOMEM,
 *         ::SEAMS_ERR_INTERNAL.
 * @par Handle access
 * Non-const because the evaluation keeps bond-weight scratch on cv. The
 * topology is not changed.
 * @par Ownership
 * The caller owns value and grad.
 * @par Thread safety
 * Exclusive use of cv. Read-only use of sys.
 * @par Example
 * @code
 * seams_cv *cv = NULL;
 * seams_cv_create(0.3, 2.8, 6, 0, 2.9, 3.0, 6, &cv);
 * seams_cv_rebuild(cv, sys, 0, NULL, 0, 0.1);
 * double n6 = 0;
 * seams_cv_value(cv, sys, SEAMS_CV_RING_COUNT, 6, &n6, NULL, 0);
 * seams_cv_destroy(cv);
 * @endcode
 */
SEAMS_API seams_status seams_cv_value(seams_cv *cv, const seams_system *sys, int kind,
                                      int arg, double *value, double *grad,
                                      size_t grad_len);

/**
 * @brief Add a soft threshold factor to a state.
 *
 * A state is the product of its factors. A factor on a collective
 * variable v is sigma((v - threshold) / width) when above is nonzero and
 * 1 - sigma((v - threshold) / width) otherwise, sigma(z) = 1 / (1 + e^-z).
 * A factor of kind ::SEAMS_CV_STATE reads an earlier state chi; with width
 * zero it is chi itself, or 1 - chi when above is zero, so states compose
 * by product (and), complement (not), and 1 - (1 - a)(1 - b) (or).
 * @param cv Workspace handle.
 * @param state Index of the state to extend, or the current number of
 *        states to open a new one.
 * @param kind A ::seams_cv_kind code.
 * @param arg Argument of kind, as for seams_cv_value(); for
 *        ::SEAMS_CV_STATE an index below state.
 * @param threshold Threshold, in the unit of the variable.
 * @param width Positive logistic width, in the unit of the variable. Zero
 *        only for ::SEAMS_CV_STATE.
 * @param above Nonzero for a factor that rises through the threshold.
 * @return ::SEAMS_OK, or ::SEAMS_ERR_NULL, ::SEAMS_ERR_RANGE,
 *         ::SEAMS_ERR_NOMEM.
 * @par Handle access
 * Appends to the state list of cv.
 * @par Ownership
 * The workspace keeps the clause until it is destroyed.
 * @par Thread safety
 * Exclusive use of cv.
 * @par Example
 * @code
 * seams_cv *cv = NULL;
 * seams_cv_create(0.3, 2.8, 6, 0, 2.9, 3.0, 6, &cv);
 * seams_cv_state_clause(cv, 0, SEAMS_CV_RING_COUNT, 6, 40.0, 4.0, 1);
 * seams_cv_destroy(cv);
 * @endcode
 */
SEAMS_API seams_status seams_cv_state_clause(seams_cv *cv, int state, int kind, int arg,
                                             double threshold, double width, int above);

/**
 * @brief Soft indicator of every state, with gradients on request.
 * @param cv Workspace handle after seams_cv_rebuild().
 * @param sys Frame with the atom count and cell of the rebuild.
 * @param chi Caller buffer of one double per state, each in [0, 1], or
 *        NULL to query n_states.
 * @param chi_len Length of chi. Ignored when chi is NULL.
 * @param n_states Receives the number of states.
 * @param grad Caller buffer of n_states rows of 3 * natoms doubles,
 *        d chi / d x in 1/Angstrom, or NULL to skip it.
 * @param grad_len Length of grad. Ignored when grad is NULL.
 * @return ::SEAMS_OK, or ::SEAMS_ERR_NULL, ::SEAMS_ERR_STATE when cv
 *         needs a rebuild for sys, ::SEAMS_ERR_RANGE, also when a chi or
 *         a gradient entry is outside double range and neither chi nor
 *         grad is written, ::SEAMS_ERR_NOMEM, ::SEAMS_ERR_INTERNAL.
 * @par Handle access
 * Non-const because the evaluation keeps bond-weight scratch on cv.
 * @par Ownership
 * The caller owns chi, n_states and grad.
 * @par Thread safety
 * Exclusive use of cv. Read-only use of sys.
 * @par Example
 * @code
 * seams_cv *cv = NULL;
 * seams_cv_create(0.3, 2.8, 6, 0, 2.9, 3.0, 6, &cv);
 * seams_cv_rebuild(cv, sys, 0, NULL, 0, 0.1);
 * seams_cv_state_clause(cv, 0, SEAMS_CV_RING_COUNT, 6, 40.0, 4.0, 1);
 * seams_cv_state_clause(cv, 1, SEAMS_CV_STATE, 0, 0.0, 0.0, 0);
 * double chi[2];
 * size_t n_states = 0;
 * seams_cv_states(cv, sys, chi, 2, &n_states, NULL, 0);
 * seams_cv_destroy(cv);
 * @endcode
 */
SEAMS_API seams_status seams_cv_states(seams_cv *cv, const seams_system *sys, double *chi,
                                       size_t chi_len, size_t *n_states, double *grad,
                                       size_t grad_len);

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
 * @param z_mid Bin centers, or NULL to skip the column.
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
 * @brief Closed cages of one ring census with an explicit face-edge window.
 *
 * Same search as seams_signature_cages(). edge NULL or "" keeps the default
 * window 2.65,3.15 in the coordinate units of the file. "none" leaves the
 * cage bound to the face census.
 *
 * @param sys Frame handle.
 * @param cutoff Positive finite cutoff, Angstrom.
 * @param type Species filter.
 * @param graph ::SEAMS_GRAPH_CUTOFF or ::SEAMS_GRAPH_MUTUAL.
 * @param spec Signature text such as "512" or "hc". Must be non-NULL.
 * @param edge Face-edge window, or NULL. "none" turns the window off.
 * @param n_cages Receives the cage count.
 * @param n_atoms Receives the vertex count.
 * @return Same codes as seams_cages_signature().
 * @par Handle access
 * Non-const only because the call clears the handle's error slot on
 * entry, so seams_handle_error() reports this call. It reads the frame
 * and keeps no scratch on the handle.
 * @par Ownership
 * The caller owns sys; the call borrows it. spec and edge are read during
 * the call only. n_cages and n_atoms are caller storage that the call fills.
 * The library keeps no pointer to any argument.
 * @par Thread safety
 * Exclusive use of sys.
 * @par Example
 * @code
 * size_t n_cages = 0, n_atoms = 0;
 * seams_signature_cages_edge(sys, 3.5, 0, SEAMS_GRAPH_CUTOFF, "512", "none", &n_cages, &n_atoms);
 * @endcode
 */
SEAMS_API seams_status seams_signature_cages_edge(seams_system *sys, double cutoff, int type,
                                                  int graph, const char *spec, const char *edge,
                                                  size_t *n_cages, size_t *n_atoms);

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
 * @param r Bin centers, or NULL to query.
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
 * @param r Bin centers, Angstrom, or NULL to query.
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

/** @brief Hydrogen-bond edge record for one donor-acceptor pair. */
typedef struct seams_hbond_edge {
  uint32_t size;     /**< Caller-set sizeof(seams_hbond_edge) on every element. */
  uint32_t donor;    /**< Donor heavy-atom index in the frame. */
  uint32_t acceptor; /**< Acceptor heavy-atom index in the frame. */
  uint32_t hydrogen; /**< Donor hydrogen index, or UINT32_MAX when omitted. */
  double d_DA;       /**< Donor-acceptor distance in Angstrom, or NaN when omitted. */
  double angle;      /**< O-O-H angle in degrees, or NaN when omitted. */
} seams_hbond_edge;

enum {
  SEAMS_HBOND_EDGE_WITH_H = 1,       /**< Populate seams_hbond_edge.hydrogen. */
  SEAMS_HBOND_EDGE_WITH_GEOMETRY = 2 /**< Populate d_DA and angle. */
};

/**
 * @brief Hydrogen-bond edges for one frame.
 *
 * The pair list is stable and deterministic: donor, then acceptor, then
 * hydrogen index. donor and acceptor are frame indices over the stored atoms.
 * When edges is non-NULL, every element must start with size set by the
 * caller to sizeof(seams_hbond_edge); the library validates size and writes
 * only up to that prefix per element.
 *
 * @param sys Frame handle.
 * @param cutoff Donor-acceptor cutoff, Angstrom.
 * @param oxygen_type Heavy-atom species filter as in seams_hbonds().
 * @param hydrogen_type Hydrogen species (positive).
 * @param dist_cutoff Acceptor-H distance, Angstrom.
 * @param angle_cutoff O-O-H angle, degrees.
 * @param fields Bitmask: SEAMS_HBOND_EDGE_WITH_H and/or
 *        SEAMS_HBOND_EDGE_WITH_GEOMETRY.
 * @param edges Caller buffer with caller-set size on each element, or NULL to query.
 * @param cap Capacity of @p edges. Ignored when @p edges is NULL.
 * @param n_out Receives the edge count.
 * @return ::SEAMS_OK, or the same error classes as seams_hbonds(), plus
 *         ::SEAMS_ERR_ARG when fields is unknown or size is invalid, and
 *         ::SEAMS_ERR_RANGE when cap is smaller than *n_out.
 * @par Handle access
 * Non-const only because the call clears the handle's error slot on
 * entry, so seams_handle_error() reports this call. It reads the frame
 * and keeps no scratch on the handle.
 * @par Ownership
 * The caller owns sys, edges, and n_out. The call borrows sys and writes
 * caller buffers only.
 * @par Thread safety
 * Exclusive use of sys.
 */
SEAMS_API seams_status seams_hbond_edges(seams_system *sys, double cutoff, int oxygen_type,
                                         int hydrogen_type, double dist_cutoff,
                                         double angle_cutoff, int fields,
                                         seams_hbond_edge *edges, size_t cap, size_t *n_out);

/** @brief One directed hydrogen bond: donor, hydrogen, acceptor. */
typedef struct seams_directed_hbond {
  uint32_t size;      /**< Caller-set sizeof(seams_directed_hbond) on every element. */
  uint32_t donor;     /**< Donor heavy-atom index in the frame. */
  uint32_t hydrogen;  /**< Hydrogen index in the frame. */
  uint32_t acceptor;  /**< Acceptor heavy-atom index in the frame. */
  double d_AH;        /**< Acceptor-hydrogen minimum-image distance, Angstrom. */
  double d_DA;        /**< Donor-acceptor minimum-image distance, Angstrom. */
  double angle;       /**< Acceptor-first angle, degrees. Linear is 0. */
} seams_directed_hbond;

/**
 * @brief Store a heavy-atom cloud and a separate hydrogen cloud.
 *
 * Heavies occupy frame indices [0, n_heavy). Hydrogens follow, in the
 * order given, so a hydrogen cloud index h is frame index n_heavy + h.
 * A null id pointer numbers that cloud from 1 across the whole frame.
 * A null mol pointer sets each molID to that atom's id, which pairs no
 * hydrogen with a heavy atom. Pass molecule ids to attach each hydrogen
 * to a heavy atom. A failed call leaves the previous frame in place.
 * @param sys Frame handle.
 * @param heavy_xyz Interleaved heavy coordinates, Angstrom.
 * @param n_heavy Heavy-atom count. Must be positive.
 * @param heavy_types One species per heavy atom.
 * @param heavy_ids One id per heavy atom, or NULL for 1..n_heavy.
 * @param heavy_mols One molID per heavy atom, or NULL to copy the atom id.
 * @param hydrogen_xyz Interleaved hydrogen coordinates, Angstrom.
 * @param n_hydrogen Hydrogen count. Zero allows the hydrogen pointers to be NULL.
 * @param hydrogen_types One species per hydrogen, required when n_hydrogen is positive.
 * @param hydrogen_ids One id per hydrogen, or NULL to continue after the heavies.
 * @param hydrogen_mols One molID per hydrogen, or NULL to copy the atom id.
 * @param box9 ox, oy, oz, lx, ly, lz, xy, xz, yz.
 * @param box_len Length of box9. Must be 9.
 * @return ::SEAMS_OK, or ::SEAMS_ERR_NULL, ::SEAMS_ERR_RANGE,
 *         ::SEAMS_ERR_BOX, ::SEAMS_ERR_ARG, ::SEAMS_ERR_NOMEM,
 *         ::SEAMS_ERR_INTERNAL.
 * @par Handle access
 * The call replaces the frame on success and clears the handle error on entry.
 * @par Ownership
 * The caller owns sys and every buffer. The call copies the values and keeps no pointer.
 * @par Thread safety
 * Exclusive use of sys.
 * @par Example
 * @code
 * const double heavy[] = {0, 0, 0, 2.8, 0, 0};
 * const double hydrogen[] = {1, 0, 0};
 * const int heavy_types[] = {8, 8};
 * const int heavy_mols[] = {1, 2};
 * const int hydrogen_types[] = {2};
 * const int hydrogen_mols[] = {1};
 * const double box9[] = {0, 0, 0, 20, 20, 20, 0, 0, 0};
 * seams_system_set_hydrogen_frame(sys, heavy, 2, heavy_types, NULL, heavy_mols,
 *                                 hydrogen, 1, hydrogen_types, NULL, hydrogen_mols,
 *                                 box9, 9);
 * @endcode
 */
SEAMS_API seams_status seams_system_set_hydrogen_frame(
    seams_system *sys, const double *heavy_xyz, size_t n_heavy, const int *heavy_types,
    const int *heavy_ids, const int *heavy_mols, const double *hydrogen_xyz,
    size_t n_hydrogen, const int *hydrogen_types, const int *hydrogen_ids,
    const int *hydrogen_mols, const double *box9, size_t box_len);

/**
 * @brief Directed hydrogen-bond edges for selected donor and acceptor types.
 *
 * Each edge is one accepting hydrogen. Rows are sorted by donor, then
 * hydrogen, then acceptor, as frame indices. The angle is acceptor-first.
 * The donor-acceptor vector is the minimum image of acceptor minus donor.
 * The acceptor-hydrogen vector is the minimum image of acceptor minus the
 * hydrogen. A linear donor-hydrogen-acceptor geometry has angle 0 degrees.
 * The acceptor-hydrogen distance is kept only when it is strictly below
 * dist_cutoff. The angle is kept when it is less than or equal to
 * angle_cutoff. The pair list keeps a donor and acceptor only when their
 * minimum image is strictly below pair_cutoff, on an orthorhombic or
 * triclinic cell. For a triclinic cell, every face separation must be at
 * least twice pair_cutoff; a thinner cell can miss valid pairs. Each
 * hydrogen is assigned to the nearest heavy atom that shares its molID,
 * then kept only when that atom's type is in donor_types. The lower index
 * wins a distance tie.
 * exclude_same_mol of 1 drops a pair that shares a molID. An empty donor
 * set or acceptor set writes zero edges. Atoms whose type equals
 * hydrogen_type are hydrogens and are never donors or acceptors.
 * Water-water and water-framework counts are two calls with different
 * type sets. When edges is non-NULL, every element must start with size
 * set to sizeof(seams_directed_hbond).
 * @param sys Frame handle.
 * @param pair_cutoff Positive finite donor-acceptor cutoff, Angstrom.
 * @param dist_cutoff Positive finite acceptor-hydrogen cutoff, Angstrom.
 * @param angle_cutoff Positive finite acceptor-first angle, degrees.
 * @param donor_types Species that may donate. NULL only when the count is 0.
 * @param n_donor_types Length of donor_types. Zero writes no edges.
 * @param acceptor_types Species that may accept. NULL only when the count is 0.
 * @param n_acceptor_types Length of acceptor_types. Zero writes no edges.
 * @param hydrogen_type Positive hydrogen species, distinct from both type sets.
 * @param exclude_same_mol 1 drops a shared molID, 0 keeps it.
 * @param edges Caller buffer, or NULL to query the count.
 * @param cap Capacity of edges. Ignored when edges is NULL.
 * @param n_out Receives the edge count.
 * @return ::SEAMS_OK, or ::SEAMS_ERR_NULL, ::SEAMS_ERR_STATE,
 *         ::SEAMS_ERR_RANGE, ::SEAMS_ERR_ARG, ::SEAMS_ERR_NOMEM,
 *         ::SEAMS_ERR_INTERNAL.
 * @par Handle access
 * Non-const only because the call clears the handle's error slot on
 * entry, so seams_handle_error() reports this call. It reads the frame
 * and keeps no scratch on the handle.
 * @par Ownership
 * The caller owns sys, the type arrays, edges, and n_out. The call borrows
 * them and writes caller buffers only.
 * @par Thread safety
 * Exclusive use of sys.
 * @par Example
 * @code
 * const int donors[] = {8};
 * const int acceptors[] = {10};
 * size_t n_out = 0;
 * seams_directed_hbonds(sys, 3.5, 2.42, 30.0, donors, 1, acceptors, 1, 2, 1,
 *                       NULL, 0, &n_out);
 * @endcode
 */
SEAMS_API seams_status seams_directed_hbonds(
    seams_system *sys, double pair_cutoff, double dist_cutoff, double angle_cutoff,
    const int *donor_types, size_t n_donor_types, const int *acceptor_types,
    size_t n_acceptor_types, int hydrogen_type, int exclude_same_mol,
    seams_directed_hbond *edges, size_t cap, size_t *n_out);

/**
 * @brief Run CHILL+, ring census, and H-bonds on raw and quenched frames.
 *
 * The callback receives one frame and returns quenched coordinates. This call
 * then runs CHILL+, primitive rings, and hydrogen-bond analysis on both raw
 * and quenched frames and emits their per-frame differences.
 * @param sys Frame handle.
 * @param callback Quench callback that writes output coordinates.
 * @param user_data Opaque callback payload.
 * @param cutoff Positive finite neighbour cutoff, Angstrom.
 * @param type Species filter for neighbour, CHILL+, and rings.
 * @param graph ::SEAMS_GRAPH_CUTOFF or ::SEAMS_GRAPH_MUTUAL.
 * @param max_depth Primitive-ring depth, from 3 through 12.
 * @param oxygen_type Heavy-atom species filter for H-bonds.
 * @param hydrogen_type Hydrogen species filter for H-bonds.
 * @param dist_cutoff Acceptor-H distance, Angstrom.
 * @param angle_cutoff O-O-H angle, degrees.
 * @param raw_labels Raw-frame CHILL+ labels, length at least atom count.
 * @param raw_labels_len Length of @p raw_labels.
 * @param quenched_labels Quenched-frame CHILL+ labels, length at least atom count.
 * @param quenched_labels_len Length of @p quenched_labels.
 * @param counts_raw Raw-frame CHILL+ population counts.
 * @param counts_quenched Quenched-frame CHILL+ population counts.
 * @param rings_raw Raw-frame primitive-ring histogram.
 * @param rings_quenched Quenched-frame primitive-ring histogram.
 * @param ring_counts_len Length of @p rings_raw and @p rings_quenched.
 * @param hbond_raw Receives raw-frame H-bond edge count.
 * @param hbond_quenched Receives quenched-frame H-bond edge count.
 * @param hbond_edge_delta Receives symmetric-difference edge count between raw
 *        and quenched H-bond edge sets.
 * @param ring_count_delta Receives sum over ring sizes of absolute raw-minus-
 *        quenched count differences.
 * @param class_change_count Receives the atom count whose CHILL+ class changed.
 * @return ::SEAMS_OK, or ::SEAMS_ERR_NULL, ::SEAMS_ERR_STATE,
 *         ::SEAMS_ERR_RANGE, ::SEAMS_ERR_BOX, ::SEAMS_ERR_ARG, ::SEAMS_ERR_NOMEM,
 *         ::SEAMS_ERR_INTERNAL.
 * @par Handle access
 * Non-const only because the call clears the handle's error slot on
 * entry, so seams_handle_error() reports this call. It reads the frame
 * and keeps no scratch on the handle.
 * @par Ownership
 * The caller owns every output buffer, callback state, and sys. The call
 * borrows callback pointers and writes caller buffers only.
 * @par Thread safety
 * Exclusive use of sys.
 */
SEAMS_API seams_status seams_quenched_analysis(
    seams_system *sys, seams_quench_callback callback, void *user_data, double cutoff, int type,
    int graph, int max_depth, int oxygen_type, int hydrogen_type, double dist_cutoff,
    double angle_cutoff, int *raw_labels, size_t raw_labels_len, int *quenched_labels,
    size_t quenched_labels_len, int counts_raw[9], int counts_quenched[9], int *rings_raw,
    int *rings_quenched, size_t ring_counts_len, int *hbond_raw, int *hbond_quenched,
    int *hbond_edge_delta, int *ring_count_delta, int *class_change_count);

/**
 * @brief Number density along one cell axis.
 *
 * z is Angstrom. rho is number density, atoms per cubic Angstrom.
 *
 * @param sys Frame handle.
 * @param type Species filter. Zero keeps every atom.
 * @param bins Bin count. Zero selects about one bin per 0.1 Angstrom.
 * @param axis 0, 1, or 2.
 * @param z Bin centers, or NULL to query.
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
 * @param color_types Non-zero colors vertices by LAMMPS type.
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
                                         int graph, int hops, int color_types,
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
 * free(key);
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
 * A lattice coefficient at or above 2^52, or a cell edge past the Euclidean
 * range of binary64, returns SEAMS_ERR_RANGE. The fractional wrap is not
 * substituted.
 *
 * @param sys Frame handle.
 * @param site_spec Text such as "1=cationHead,2=anion". Must be non-NULL.
 * @param n_pairs Receives the contact count, or NULL.
 * @param n_cation Receives the cation molecule count, or NULL.
 * @param n_anion Receives the anion molecule count, or NULL. At least
 *        one of the three outputs must be non-NULL.
 * @return ::SEAMS_OK, or ::SEAMS_ERR_NULL, ::SEAMS_ERR_STATE,
 *         ::SEAMS_ERR_RANGE, ::SEAMS_ERR_ARG, ::SEAMS_ERR_NOMEM,
 *         ::SEAMS_ERR_INTERNAL.
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
 * if (seams_system_natoms(sys, &n) == SEAMS_OK) {
 *   labels = malloc(n * sizeof *labels);
 *   free(labels);
 * }
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
 * free(sym);
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
 * free(xyz);
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

/**
 * @defgroup capi_traj Trajectories and dynamics
 * @brief Periodic unwrapping, molecular dipoles, and streaming time correlations.
 * @{
 */

/**
 * @brief Copy the stored molecule ids out.
 * @param sys Frame handle.
 * @param mols Caller buffer for n ints, in stored order.
 * @param len Length of mols. Must be at least n.
 * @return ::SEAMS_OK; ::SEAMS_ERR_NULL; ::SEAMS_ERR_STATE when no frame
 *         is stored; ::SEAMS_ERR_RANGE for a short buffer.
 * @par Ownership
 * The caller owns mols.
 * @par Thread safety
 * Read-only: may run alongside other const calls on sys.
 * @par Example
 * @code
 * int mols[8] = {0};
 * seams_mols(sys, mols, 8);
 * @endcode
 */
SEAMS_API seams_status seams_mols(const seams_system *sys, int *mols, size_t len);

/**
 * @brief Replace the stored molecule ids on a loaded frame.
 * @param sys Frame handle.
 * @param mols Caller buffer of n molecule ids.
 * @param len Length of mols. Must equal the stored atom count.
 * @return ::SEAMS_OK; ::SEAMS_ERR_NULL; ::SEAMS_ERR_STATE when no frame
 *         is stored; ::SEAMS_ERR_RANGE when len differs from the atom count.
 * @par Handle access
 * Non-const: updates stored molecule ids and clears the handle error slot.
 * @par Ownership
 * The library copies mols and keeps no pointer to caller memory.
 * @par Thread safety
 * Exclusive use of sys.
 * @par Example
 * @code
 * const int mols[8] = {1, 1, 2, 2, 3, 3, 4, 4};
 * seams_set_mols(sys, mols, 8);
 * @endcode
 */
SEAMS_API seams_status seams_set_mols(seams_system *sys, const int *mols, size_t len);

/**
 * @brief Intact periodic centres of mass or geometry for molecules in a frame.
 *
 * Atoms sharing a molecule id are assembled onto the first atom of that
 * molecule through the minimum image of the stored orthorhombic or triclinic
 * cell before averaging. Pass NULL for both mol_ids_out and centres_out to
 * query the molecule count in n_mols.
 *
 * @param sys Frame handle.
 * @param mols Optional molecule id per atom, or NULL to use stored molecule ids.
 * @param mols_len Length of mols when mols is non-NULL; must equal n, or 0 when NULL.
 * @param masses Optional positive weight per atom, or NULL for equal weights.
 * @param masses_len Length of masses when masses is non-NULL; must equal n, or 0 when NULL.
 * @param mol_ids_out Receives one molecule id per molecule, or NULL.
 * @param centres_out Receives 3 * n_mols coordinates in Angstrom, or NULL.
 * @param cap Molecule capacity of non-NULL output buffers.
 * @param n_mols Receives the molecule count.
 * @return ::SEAMS_OK; ::SEAMS_ERR_NULL; ::SEAMS_ERR_STATE; ::SEAMS_ERR_RANGE;
 *         ::SEAMS_ERR_ARG; ::SEAMS_ERR_BOX; ::SEAMS_ERR_NOMEM.
 * @par Ownership
 * The caller owns all buffers.
 * @par Thread safety
 * Read-only on sys.
 * @par Example
 * @code
 * size_t nm = 0;
 * seams_molecule_centres(sys, NULL, 0, NULL, 0, NULL, NULL, 0, &nm);
 * @endcode
 */
SEAMS_API seams_status seams_molecule_centres(const seams_system *sys,
                                              const int *mols, size_t mols_len,
                                              const double *masses,
                                              size_t masses_len,
                                              int *mol_ids_out,
                                              double *centres_out, size_t cap,
                                              size_t *n_mols);

/**
 * @brief Unwrap periodic coordinates across consecutive frames.
 *
 * Each step adds the minimum-image displacement under an orthorhombic or
 * triclinic nine-number cell (ox, oy, oz, lx, ly, lz, xy, xz, yz). Pass
 * nboxes equal to 1 for a constant cell or nframes for a deforming cell.
 *
 * @param xyz Wrapped coordinates, nframes * 3 * natoms doubles in Angstrom.
 * @param nframes Number of frames, at least 1.
 * @param natoms Number of atoms per frame, at least 1.
 * @param boxes Cell array of 9 * nboxes doubles.
 * @param nboxes Cell count, either 1 or nframes.
 * @param unwrapped Output buffer of unwrapped_len doubles.
 * @param unwrapped_len Capacity in doubles, at least nframes * 3 * natoms.
 * @return ::SEAMS_OK; ::SEAMS_ERR_NULL; ::SEAMS_ERR_RANGE; ::SEAMS_ERR_BOX;
 *         ::SEAMS_ERR_ARG; ::SEAMS_ERR_NOMEM.
 * @par Ownership
 * The caller owns xyz, boxes, and unwrapped.
 * @par Thread safety
 * Pure aside from the thread-local error string.
 * @par Example
 * @code
 * const double box[9] = {0, 0, 0, 10, 10, 10, 0, 0, 0};
 * const double wrapped[6] = {9.5, 1.0, 1.0, 0.5, 1.0, 1.0};
 * double unwrapped[6] = {0};
 * seams_traj_unwrap_positions(wrapped, 2, 1, box, 1, unwrapped, 6);
 * @endcode
 */
SEAMS_API seams_status seams_traj_unwrap_positions(const double *xyz,
                                                   size_t nframes,
                                                   size_t natoms,
                                                   const double *boxes,
                                                   size_t nboxes,
                                                   double *unwrapped,
                                                   size_t unwrapped_len);

/**
 * @brief Assemble per-molecule centres and unwrap them across frames.
 *
 * Pass NULL for both mol_ids_out and centres_out to query n_mols first.
 *
 * @param xyz Wrapped atom coordinates, nframes * 3 * natoms doubles.
 * @param nframes Number of frames, at least 1.
 * @param natoms Number of atoms per frame, at least 1.
 * @param boxes Cell array of 9 * nboxes doubles.
 * @param nboxes Cell count, either 1 or nframes.
 * @param mols Molecule id per atom, length mols_len.
 * @param mols_len Must equal natoms.
 * @param masses Optional positive mass per atom, or NULL for equal weights.
 * @param masses_len Must equal natoms when masses is non-NULL, or 0 when NULL.
 * @param mol_ids_out Receives unique molecule ids, or NULL.
 * @param mol_ids_cap Capacity of mol_ids_out in ints when non-NULL.
 * @param centres_out Receives nframes * 3 * n_mols unwrapped coordinates, or NULL.
 * @param centres_len Capacity of centres_out in doubles when non-NULL.
 * @param n_mols Receives the molecule count.
 * @return ::SEAMS_OK; ::SEAMS_ERR_NULL; ::SEAMS_ERR_RANGE; ::SEAMS_ERR_BOX;
 *         ::SEAMS_ERR_ARG; ::SEAMS_ERR_NOMEM.
 * @par Ownership
 * The caller owns all buffers.
 * @par Thread safety
 * Pure aside from the thread-local error string.
 * @par Example
 * @code
 * const double box[9] = {0, 0, 0, 10, 10, 10, 0, 0, 0};
 * const double xyz[12] = {9.8, 0, 0, 0.2, 0, 0, 0.4, 0, 0, 0.8, 0, 0};
 * const int mols[2] = {1, 1};
 * size_t nm = 0;
 * double out[6] = {0};
 * seams_traj_unwrap_centres(xyz, 2, 2, box, 1, mols, 2, NULL, 0, NULL, 0, out, 6, &nm);
 * @endcode
 */
SEAMS_API seams_status seams_traj_unwrap_centres(
    const double *xyz, size_t nframes, size_t natoms, const double *boxes,
    size_t nboxes, const int *mols, size_t mols_len, const double *masses,
    size_t masses_len, int *mol_ids_out, size_t mol_ids_cap,
    double *centres_out, size_t centres_len, size_t *n_mols);

/**
 * @brief Water orientation vectors: O to H-midpoint and H-to-H under PBC.
 *
 * Each water pairs one oxygen_type atom with two hydrogen_type atoms via
 * shared molecule ids when available, or via minimum-image covalent distance
 * otherwise. Pass NULL for all five output arrays to query n_waters.
 *
 * @param sys Frame handle.
 * @param oxygen_type Species integer of water oxygen.
 * @param hydrogen_type Species integer of water hydrogen, distinct from oxygen_type.
 * @param mol_ids Receives molecule ids, length cap, or NULL.
 * @param oh_mid Receives O-to-H-midpoint vectors in Angstrom, 3 * cap doubles, or NULL.
 * @param bisector Receives unit bisector vectors, 3 * cap doubles, or NULL.
 * @param hh Receives H1-to-H2 vectors in Angstrom, 3 * cap doubles, or NULL.
 * @param hh_unit Receives unit H1-to-H2 vectors, 3 * cap doubles, or NULL.
 * @param cap Molecule capacity of every non-NULL output array.
 * @param n_waters Receives the water count.
 * @return ::SEAMS_OK; ::SEAMS_ERR_NULL; ::SEAMS_ERR_STATE; ::SEAMS_ERR_RANGE;
 *         ::SEAMS_ERR_ARG; ::SEAMS_ERR_BOX; ::SEAMS_ERR_NOMEM.
 * @par Ownership
 * The caller owns all output arrays.
 * @par Thread safety
 * Read-only on sys.
 * @par Example
 * @code
 * size_t nw = 0;
 * seams_water_vectors(sys, 1, 2, NULL, NULL, NULL, NULL, NULL, 0, &nw);
 * @endcode
 */
SEAMS_API seams_status seams_water_vectors(const seams_system *sys,
                                           int oxygen_type, int hydrogen_type,
                                           int *mol_ids, double *oh_mid,
                                           double *bisector, double *hh,
                                           double *hh_unit, size_t cap,
                                           size_t *n_waters);

/**
 * @brief Generic per-molecule dipole or geometric orientation vectors.
 *
 * With charges, computes sum_a q_a (r_a - r_centre) in elementary-charge
 * Angstrom on each PBC-assembled molecule. Without charges, computes the mean
 * minimum-image vector from each molecule's first atom to its remaining atoms.
 *
 * @param sys Frame handle.
 * @param mols Optional molecule id per atom, or NULL to use stored molecule ids.
 * @param mols_len Length of mols when non-NULL, or 0 when NULL.
 * @param charges Optional partial charge per atom, or NULL for geometric mode.
 * @param charges_len Length of charges when non-NULL, or 0 when NULL.
 * @param mol_ids_out Receives molecule ids, length cap, or NULL.
 * @param dipoles Receives 3 * cap raw dipole vectors, or NULL.
 * @param unit_vectors Receives 3 * cap normalized orientation vectors, or NULL.
 * @param cap Molecule capacity of non-NULL output arrays.
 * @param n_mols Receives the molecule count.
 * @return ::SEAMS_OK; ::SEAMS_ERR_NULL; ::SEAMS_ERR_STATE; ::SEAMS_ERR_RANGE;
 *         ::SEAMS_ERR_ARG; ::SEAMS_ERR_BOX; ::SEAMS_ERR_NOMEM.
 * @par Ownership
 * The caller owns all buffers.
 * @par Thread safety
 * Read-only on sys.
 * @par Example
 * @code
 * size_t nm = 0;
 * seams_molecule_dipoles(sys, NULL, 0, NULL, 0, NULL, NULL, NULL, 0, &nm);
 * @endcode
 */
SEAMS_API seams_status seams_molecule_dipoles(
    const seams_system *sys, const int *mols, size_t mols_len,
    const double *charges, size_t charges_len, int *mol_ids_out,
    double *dipoles, double *unit_vectors, size_t cap, size_t *n_mols);

/**
 * @brief Project molecular vectors onto an axis for P1 and P2 order parameters.
 *
 * Computes cos_theta = u . e and P2 = 0.5 * (3 cos_theta^2 - 1) for each
 * non-zero vector, where e is the normalized axis. Zero-length vectors receive
 * NaN in per-molecule outputs and are omitted from mean_p1 and mean_p2.
 *
 * @param vectors Input 3-vectors, 3 * n_vectors doubles (or NULL when n_vectors is 0).
 * @param n_vectors Number of 3-vectors.
 * @param axis Non-zero reference axis of 3 finite doubles.
 * @param cos_theta Optional output buffer of cap doubles for cos(theta), or NULL.
 * @param p2 Optional output buffer of cap doubles for P2(cos(theta)), or NULL.
 * @param cap Capacity of non-NULL per-vector buffers; must be at least n_vectors.
 * @param mean_p1 Receives ensemble mean of cos(theta), or NULL.
 * @param mean_p2 Receives ensemble mean of P2(cos(theta)), or NULL.
 * @return ::SEAMS_OK; ::SEAMS_ERR_NULL; ::SEAMS_ERR_RANGE; ::SEAMS_ERR_ARG;
 *         ::SEAMS_ERR_NOMEM.
 * @par Ownership
 * The caller owns all arrays.
 * @par Thread safety
 * Pure aside from the thread-local error string.
 * @par Example
 * @code
 * const double vecs[6] = {0, 0, 1, 0, 0, 2};
 * const double axis[3] = {0, 0, 1};
 * double p1 = 0, p2 = 0;
 * seams_orientation_order(vecs, 2, axis, NULL, NULL, 0, &p1, &p2);
 * @endcode
 */
SEAMS_API seams_status seams_orientation_order(const double *vectors,
                                               size_t n_vectors,
                                               const double axis[3],
                                               double *cos_theta, double *p2,
                                               size_t cap, double *mean_p1,
                                               double *mean_p2);

/**
 * @brief Create a streaming multi-origin MSD accumulator with axial/lateral split.
 *
 * @param natoms Particle or molecule count per frame, at least 1. Above
 *        PTRDIFF_MAX / (3 * sizeof(double)) the coordinates cannot be
 *        allocated, which is ::SEAMS_ERR_RANGE.
 * @param max_lag Maximum lag in frames, or 0 to retain all lags. Above
 *        PTRDIFF_MAX / sizeof(double) - 1 the lag bins cannot be allocated,
 *        which is ::SEAMS_ERR_RANGE.
 * @param origin_stride Frame stride between time origins, at least 1.
 * @param dt Positive finite frame spacing in physical time units.
 * @param axis Non-zero 3-vector defining the axial direction.
 * @param out Receives the new handle; set to NULL on failure.
 * @return ::SEAMS_OK; ::SEAMS_ERR_NULL; ::SEAMS_ERR_RANGE; ::SEAMS_ERR_ARG;
 *         ::SEAMS_ERR_NOMEM.
 * @par Ownership
 * The caller frees *out with seams_msd_destroy().
 * @par Thread safety
 * Pure aside from handle registration and the thread-local error string.
 * @par Example
 * @code
 * const double axis[3] = {0, 0, 1};
 * seams_msd *msd = NULL;
 * seams_msd_create(4, 16, 1, 0.1, axis, &msd);
 * seams_msd_destroy(msd);
 * @endcode
 */
SEAMS_API seams_status seams_msd_create(size_t natoms, size_t max_lag,
                                        size_t origin_stride, double dt,
                                        const double axis[3], seams_msd **out);

/**
 * @brief Free an MSD accumulator.
 * @param msd Handle from seams_msd_create(), or NULL.
 * @par Ownership
 * Releases the handle and its ring buffer.
 * @par Thread safety
 * Exclusive use of msd.
 * @par Example
 * @code
 * seams_msd_destroy(NULL);
 * @endcode
 */
SEAMS_API void seams_msd_destroy(seams_msd *msd);

/**
 * @brief Push one frame into a streaming MSD accumulator.
 *
 * Pass box9 as NULL (with box_len 0) when xyz is already unwrapped, or pass a
 * nine-number cell (with box_len 9) to unwrap incrementally across frames.
 *
 * @param msd MSD handle.
 * @param xyz Coordinates of 3 * natoms doubles in Angstrom.
 * @param natoms Particle count; must match the count passed to seams_msd_create().
 * @param box9 Optional 9-number periodic cell, or NULL for unwrapped input.
 * @param box_len Must be 9 when box9 is non-NULL, or 0 when box9 is NULL.
 * @return ::SEAMS_OK; ::SEAMS_ERR_NULL; ::SEAMS_ERR_RANGE; ::SEAMS_ERR_BOX;
 *         ::SEAMS_ERR_ARG; ::SEAMS_ERR_STATE; ::SEAMS_ERR_NOMEM.
 * @par Ownership
 * The library copies coordinates into active origin slots as needed.
 * @par Thread safety
 * Exclusive use of msd.
 * @par Example
 * @code
 * const double axis[3] = {0, 0, 1};
 * const double f0[3] = {0, 0, 0};
 * seams_msd *msd = NULL;
 * seams_msd_create(1, 8, 1, 1.0, axis, &msd);
 * seams_msd_push(msd, f0, 1, NULL, 0);
 * seams_msd_destroy(msd);
 * @endcode
 */
SEAMS_API seams_status seams_msd_push(seams_msd *msd, const double *xyz,
                                      size_t natoms, const double *box9,
                                      size_t box_len);

/**
 * @brief Read accumulated lag times and total, axial, and lateral MSD curves.
 *
 * Any of lag_time, total, axial, lateral, and n_origins may be NULL to skip
 * that column; pass all five as NULL to query n_lags.
 *
 * @param msd MSD handle.
 * @param lag_time Receives lag times k * dt, or NULL.
 * @param total Receives 3D total MSD, or NULL.
 * @param axial Receives 1D axial MSD along axis, or NULL.
 * @param lateral Receives 2D lateral MSD perpendicular to axis, or NULL.
 * @param n_origins Receives time-origin counts per lag, or NULL.
 * @param cap Capacity of non-NULL output arrays.
 * @param n_lags Receives the available lag count.
 * @return ::SEAMS_OK; ::SEAMS_ERR_NULL; ::SEAMS_ERR_STATE; ::SEAMS_ERR_RANGE;
 *         ::SEAMS_ERR_NOMEM.
 * @par Ownership
 * The caller owns all output buffers.
 * @par Thread safety
 * Read-only on msd.
 * @par Example
 * @code
 * const double axis[3] = {0, 0, 1};
 * const double f0[3] = {0, 0, 0}, f1[3] = {0, 0, 1};
 * seams_msd *msd = NULL;
 * size_t nlags = 0;
 * seams_msd_create(1, 4, 1, 1.0, axis, &msd);
 * seams_msd_push(msd, f0, 1, NULL, 0);
 * seams_msd_push(msd, f1, 1, NULL, 0);
 * seams_msd_curve(msd, NULL, NULL, NULL, NULL, NULL, 0, &nlags);
 * seams_msd_destroy(msd);
 * @endcode
 */
SEAMS_API seams_status seams_msd_curve(const seams_msd *msd, double *lag_time,
                                       double *total, double *axial,
                                       double *lateral, size_t *n_origins,
                                       size_t cap, size_t *n_lags);

/**
 * @brief Linear Einstein diffusion fit on an MSD accumulator.
 *
 * Fits MSD_total = 6 D_total t + c, MSD_axial = 2 D_axial t + c, and
 * MSD_lateral = 4 D_lateral t + c over [t_start, t_end]. Pass t_start = 0 and
 * t_end = 0 for the default [0.1 t_max, 0.5 t_max] window.
 *
 * @param msd MSD handle with at least 2 lag points.
 * @param t_start Window start time, >= 0.
 * @param t_end Window end time, > t_start (or both 0 for default window).
 * @param d_total Receives 3D diffusion coefficient, or NULL.
 * @param d_axial Receives 1D axial diffusion coefficient, or NULL.
 * @param d_lateral Receives 2D lateral diffusion coefficient, or NULL.
 * @param alpha Receives the log-log anomalous exponent, or NULL. Receives NaN
 *        when the window holds fewer than two lags with t > 0 and
 *        MSD_total > 1e-20, so the log-log fit was skipped; test it with
 *        isnan(). The diffusion coefficients are still the linear fit.
 * @return ::SEAMS_OK when the linear fit succeeded, including a skipped
 *         log-log fit (alpha NaN); ::SEAMS_ERR_NULL; ::SEAMS_ERR_STATE;
 *         ::SEAMS_ERR_RANGE; ::SEAMS_ERR_ARG; ::SEAMS_ERR_NOMEM.
 * @par Ownership
 * Writes scalar outputs owned by the caller.
 * @par Thread safety
 * Read-only on msd.
 * @par Example
 * @code
 * const double axis[3] = {0, 0, 1};
 * const double f0[3] = {0, 0, 0}, f1[3] = {0, 0, 1}, f2[3] = {0, 0, 2};
 * seams_msd *msd = NULL;
 * double dtot = 0, dax = 0, dlat = 0, alpha = 0;
 * seams_msd_create(1, 4, 1, 1.0, axis, &msd);
 * seams_msd_push(msd, f0, 1, NULL, 0);
 * seams_msd_push(msd, f1, 1, NULL, 0);
 * seams_msd_push(msd, f2, 1, NULL, 0);
 * seams_msd_fit(msd, 0.0, 0.0, &dtot, &dax, &dlat, &alpha);
 * seams_msd_destroy(msd);
 * @endcode
 */
SEAMS_API seams_status seams_msd_fit(const seams_msd *msd, double t_start,
                                     double t_end, double *d_total,
                                     double *d_axial, double *d_lateral,
                                     double *alpha);

/**
 * @brief Create a streaming H-bond lifetime accumulator from per-frame edges.
 *
 * @param max_lag Maximum lag in frames, or 0 to retain all lags. Above
 *        PTRDIFF_MAX / sizeof(double) - 1 the lag bins cannot be allocated,
 *        which is ::SEAMS_ERR_RANGE.
 * @param origin_stride Frame stride between time origins, at least 1.
 * @param dt Positive finite frame spacing.
 * @param directed Non-zero treats (donor, acceptor) as directed; zero canonicalizes pairs.
 * @param out Receives the new handle; set to NULL on failure.
 * @return ::SEAMS_OK; ::SEAMS_ERR_NULL; ::SEAMS_ERR_RANGE; ::SEAMS_ERR_ARG;
 *         ::SEAMS_ERR_NOMEM.
 * @par Ownership
 * The caller frees *out with seams_hbond_tcf_destroy().
 * @par Thread safety
 * Pure aside from handle registration and the thread-local error string.
 * @par Example
 * @code
 * seams_hbond_tcf *tcf = NULL;
 * seams_hbond_tcf_create(16, 1, 0.1, 1, &tcf);
 * seams_hbond_tcf_destroy(tcf);
 * @endcode
 */
SEAMS_API seams_status seams_hbond_tcf_create(size_t max_lag,
                                              size_t origin_stride, double dt,
                                              int directed,
                                              seams_hbond_tcf **out);

/**
 * @brief Free an H-bond lifetime accumulator.
 * @param tcf Handle from seams_hbond_tcf_create(), or NULL.
 * @par Ownership
 * Releases the handle and its active origin sets.
 * @par Thread safety
 * Exclusive use of tcf.
 * @par Example
 * @code
 * seams_hbond_tcf_destroy(NULL);
 * @endcode
 */
SEAMS_API void seams_hbond_tcf_destroy(seams_hbond_tcf *tcf);

/**
 * @brief Push one frame's hydrogen-bond edge set into the lifetime accumulator.
 *
 * Duplicate edges within the same frame are deduplicated before updating
 * intermittent C_HB(t) and continuous survival S_HB(t).
 *
 * @param tcf H-bond TCF handle.
 * @param donors Donor ids of length n_edges (may be NULL when n_edges is 0).
 * @param acceptors Acceptor ids of length n_edges (may be NULL when n_edges is 0).
 * @param n_edges Number of edges in this frame.
 * @return ::SEAMS_OK; ::SEAMS_ERR_NULL; ::SEAMS_ERR_NOMEM.
 * @par Ownership
 * The library copies edge ids and keeps no caller pointer.
 * @par Thread safety
 * Exclusive use of tcf.
 * @par Example
 * @code
 * const int don[1] = {1}, acc[1] = {2};
 * seams_hbond_tcf *tcf = NULL;
 * seams_hbond_tcf_create(8, 1, 1.0, 1, &tcf);
 * seams_hbond_tcf_push(tcf, don, acc, 1);
 * seams_hbond_tcf_destroy(tcf);
 * @endcode
 */
SEAMS_API seams_status seams_hbond_tcf_push(seams_hbond_tcf *tcf,
                                            const int *donors,
                                            const int *acceptors,
                                            size_t n_edges);

/**
 * @brief Read intermittent C_HB(t) and continuous S_HB(t) lifetime curves.
 *
 * @param tcf H-bond TCF handle.
 * @param lag_time Receives lag times k * dt, or NULL.
 * @param c_intermittent Receives intermittent autocorrelation C_HB, or NULL.
 * @param s_continuous Receives unbroken survival probability S_HB, or NULL.
 * @param n_origins Receives origin counts per lag, or NULL.
 * @param cap Capacity of non-NULL output arrays.
 * @param n_lags Receives the available lag count.
 * @return ::SEAMS_OK; ::SEAMS_ERR_NULL; ::SEAMS_ERR_STATE; ::SEAMS_ERR_RANGE;
 *         ::SEAMS_ERR_NOMEM.
 * @par Ownership
 * The caller owns all output buffers.
 * @par Thread safety
 * Read-only on tcf.
 * @par Example
 * @code
 * const int don[1] = {1}, acc[1] = {2};
 * seams_hbond_tcf *tcf = NULL;
 * size_t nlags = 0;
 * seams_hbond_tcf_create(4, 1, 1.0, 1, &tcf);
 * seams_hbond_tcf_push(tcf, don, acc, 1);
 * seams_hbond_tcf_push(tcf, don, acc, 1);
 * seams_hbond_tcf_curve(tcf, NULL, NULL, NULL, NULL, 0, &nlags);
 * seams_hbond_tcf_destroy(tcf);
 * @endcode
 */
SEAMS_API seams_status seams_hbond_tcf_curve(const seams_hbond_tcf *tcf,
                                             double *lag_time,
                                             double *c_intermittent,
                                             double *s_continuous,
                                             size_t *n_origins, size_t cap,
                                             size_t *n_lags);

/**
 * @brief Integrated and single-exponential lifetimes for C_HB and S_HB.
 *
 * Pass t_start = 0 and t_end = 0 to fit across the full positive decay range.
 *
 * @param tcf H-bond TCF handle with at least 2 lag points.
 * @param t_start Fit window start time, >= 0.
 * @param t_end Fit window end time, > t_start (or both 0 for full range).
 * @param tau_c_int Receives trapezoidal integral of C_HB, or NULL.
 * @param tau_c_exp Receives exponential decay time of C_HB, or NULL.
 * @param tau_s_int Receives trapezoidal integral of S_HB, or NULL.
 * @param tau_s_exp Receives exponential decay time of S_HB, or NULL.
 * @return ::SEAMS_OK; ::SEAMS_ERR_NULL; ::SEAMS_ERR_STATE; ::SEAMS_ERR_RANGE;
 *         ::SEAMS_ERR_ARG; ::SEAMS_ERR_NOMEM.
 * @par Ownership
 * Writes scalar outputs owned by the caller.
 * @par Thread safety
 * Read-only on tcf.
 * @par Example
 * @code
 * const int don[2] = {1, 2}, acc[2] = {3, 4};
 * seams_hbond_tcf *tcf = NULL;
 * double tc = 0, ts = 0;
 * seams_hbond_tcf_create(4, 1, 1.0, 1, &tcf);
 * seams_hbond_tcf_push(tcf, don, acc, 2);
 * seams_hbond_tcf_push(tcf, don, acc, 1);
 * seams_hbond_tcf_fit(tcf, 0.0, 0.0, &tc, NULL, &ts, NULL);
 * seams_hbond_tcf_destroy(tcf);
 * @endcode
 */
SEAMS_API seams_status seams_hbond_tcf_fit(const seams_hbond_tcf *tcf,
                                           double t_start, double t_end,
                                           double *tau_c_int, double *tau_c_exp,
                                           double *tau_s_int,
                                           double *tau_s_exp);

/**
 * @brief Create a streaming orientational correlation accumulator for C1 and C2.
 *
 * @param n_vectors Number of molecular vectors per frame, at least 1. Above
 *        PTRDIFF_MAX / (3 * sizeof(double)) the vectors cannot be allocated,
 *        which is ::SEAMS_ERR_RANGE.
 * @param max_lag Maximum lag in frames, or 0 to retain all lags. Above
 *        PTRDIFF_MAX / sizeof(double) - 1 the lag bins cannot be allocated,
 *        which is ::SEAMS_ERR_RANGE.
 * @param origin_stride Frame stride between time origins, at least 1.
 * @param dt Positive finite frame spacing.
 * @param out Receives the new handle; set to NULL on failure.
 * @return ::SEAMS_OK; ::SEAMS_ERR_NULL; ::SEAMS_ERR_RANGE; ::SEAMS_ERR_ARG;
 *         ::SEAMS_ERR_NOMEM.
 * @par Ownership
 * The caller frees *out with seams_reorient_destroy().
 * @par Thread safety
 * Pure aside from handle registration and the thread-local error string.
 * @par Example
 * @code
 * seams_reorient *rot = NULL;
 * seams_reorient_create(4, 16, 1, 0.1, &rot);
 * seams_reorient_destroy(rot);
 * @endcode
 */
SEAMS_API seams_status seams_reorient_create(size_t n_vectors, size_t max_lag,
                                             size_t origin_stride, double dt,
                                             seams_reorient **out);

/**
 * @brief Free a reorientation accumulator.
 * @param reorient Handle from seams_reorient_create(), or NULL.
 * @par Ownership
 * Releases the handle and its ring buffer.
 * @par Thread safety
 * Exclusive use of reorient.
 * @par Example
 * @code
 * seams_reorient_destroy(NULL);
 * @endcode
 */
SEAMS_API void seams_reorient_destroy(seams_reorient *reorient);

/**
 * @brief Push one frame of molecular 3-vectors into a reorientation accumulator.
 *
 * Non-zero vectors are normalized per molecule; zero-norm vectors are skipped
 * in correlation pairs.
 *
 * @param reorient Reorientation handle.
 * @param vectors Array of 3 * n_vectors finite doubles.
 * @param n_vectors Vector count; must match seams_reorient_create().
 * @return ::SEAMS_OK; ::SEAMS_ERR_NULL; ::SEAMS_ERR_RANGE; ::SEAMS_ERR_ARG;
 *         ::SEAMS_ERR_NOMEM.
 * @par Ownership
 * The library copies normalized vectors into active origin slots.
 * @par Thread safety
 * Exclusive use of reorient.
 * @par Example
 * @code
 * const double v0[3] = {0, 0, 1};
 * seams_reorient *rot = NULL;
 * seams_reorient_create(1, 4, 1, 1.0, &rot);
 * seams_reorient_push(rot, v0, 1);
 * seams_reorient_destroy(rot);
 * @endcode
 */
SEAMS_API seams_status seams_reorient_push(seams_reorient *reorient,
                                           const double *vectors,
                                           size_t n_vectors);

/**
 * @brief Read accumulated C1(t) = <P1(u(0).u(t))> and C2(t) = <P2(u(0).u(t))>.
 *
 * @param reorient Reorientation handle.
 * @param lag_time Receives lag times k * dt, or NULL.
 * @param c1 Receives first Legendre correlation C1, or NULL.
 * @param c2 Receives second Legendre correlation C2, or NULL.
 * @param n_origins Receives origin counts per lag, or NULL.
 * @param cap Capacity of non-NULL output arrays.
 * @param n_lags Receives the available lag count.
 * @return ::SEAMS_OK; ::SEAMS_ERR_NULL; ::SEAMS_ERR_STATE; ::SEAMS_ERR_RANGE;
 *         ::SEAMS_ERR_NOMEM.
 * @par Ownership
 * The caller owns all output buffers.
 * @par Thread safety
 * Read-only on reorient.
 * @par Example
 * @code
 * const double v0[3] = {0, 0, 1}, v1[3] = {0, 1, 0};
 * seams_reorient *rot = NULL;
 * size_t nlags = 0;
 * seams_reorient_create(1, 4, 1, 1.0, &rot);
 * seams_reorient_push(rot, v0, 1);
 * seams_reorient_push(rot, v1, 1);
 * seams_reorient_curve(rot, NULL, NULL, NULL, NULL, 0, &nlags);
 * seams_reorient_destroy(rot);
 * @endcode
 */
SEAMS_API seams_status seams_reorient_curve(const seams_reorient *reorient,
                                            double *lag_time, double *c1,
                                            double *c2, size_t *n_origins,
                                            size_t cap, size_t *n_lags);

/**
 * @brief Integrated and single-exponential reorientation times tau_1 and tau_2.
 *
 * @param reorient Reorientation handle with at least 2 lag points.
 * @param t_start Fit window start time, >= 0.
 * @param t_end Fit window end time, > t_start (or both 0 for full range).
 * @param tau1_int Receives trapezoidal integral of C1, or NULL.
 * @param tau1_exp Receives exponential decay time of C1, or NULL.
 * @param tau2_int Receives trapezoidal integral of C2, or NULL.
 * @param tau2_exp Receives exponential decay time of C2, or NULL.
 * @return ::SEAMS_OK; ::SEAMS_ERR_NULL; ::SEAMS_ERR_STATE; ::SEAMS_ERR_RANGE;
 *         ::SEAMS_ERR_ARG; ::SEAMS_ERR_NOMEM.
 * @par Ownership
 * Writes scalar outputs owned by the caller.
 * @par Thread safety
 * Read-only on reorient.
 * @par Example
 * @code
 * const double v0[3] = {0, 0, 1}, v1[3] = {0.6, 0, 0.8};
 * seams_reorient *rot = NULL;
 * double t1 = 0, t2 = 0;
 * seams_reorient_create(1, 4, 1, 1.0, &rot);
 * seams_reorient_push(rot, v0, 1);
 * seams_reorient_push(rot, v1, 1);
 * seams_reorient_fit(rot, 0.0, 0.0, &t1, NULL, &t2, NULL);
 * seams_reorient_destroy(rot);
 * @endcode
 */
SEAMS_API seams_status seams_reorient_fit(const seams_reorient *reorient,
                                          double t_start, double t_end,
                                          double *tau1_int, double *tau1_exp,
                                          double *tau2_int, double *tau2_exp);

/** @} */

/**
 * @defgroup capi_pore Pore geometry
 * @brief Fractional density, cylindrical pore profiles, and host-guest contacts.
 *
 * Slab and shell volumes use the triclinic cell, det(H). A fractional
 * bin along a cell edge has volume V/n. An arbitrary direction bins the
 * projection of the reduced cell, and each bin volume is the measure of
 * that slice. Number density is atoms per cubic Angstrom. Mass density
 * is the caller's mass unit per cubic Angstrom, or mass 1 when masses
 * is NULL. A mass given in amu is therefore amu per cubic Angstrom, and
 * multiplying by 1.66054 converts that density to grams per cubic
 * centimetre. The library does not apply the factor.
 * @{
 */

/**
 * @brief Number and mass density along a cell axis or a Cartesian direction.
 *
 * coordinate is the bin centre in the shifted fractional coordinate, on
 * [0, 1). Bin 0 starts at origin. origin is wrapped into [0, 1). A
 * non-NULL direction is three Cartesian components and replaces axis.
 * Pass NULL output columns to query the bin count.
 *
 * @param sys Frame handle.
 * @param type Species filter. Zero keeps every atom.
 * @param bins Positive bin count.
 * @param axis 0, 1, or 2. Ignored when direction is not NULL.
 * @param direction Three Cartesian components, or NULL for a cell axis.
 * @param direction_len 3 when direction is not NULL. Ignored otherwise.
 * @param origin Fractional origin. Any finite value is wrapped into [0, 1).
 * @param masses One mass per stored atom, or NULL for mass 1.
 * @param masses_len Length of masses. Must equal the atom count when
 *        masses is not NULL.
 * @param coordinate Bin centres, or NULL to skip.
 * @param number_density Number density, or NULL to skip.
 * @param mass_density Mass density in the caller's mass units per cubic
 *        Angstrom, or NULL to skip. An amu mass gives amu per cubic
 *        Angstrom; times 1.66054 that is grams per cubic centimetre.
 * @param bin_volume Bin volume in cubic Angstrom, or NULL to skip.
 * @param cap Capacity of each non-NULL column.
 * @param n_out Receives the bin count.
 * @return ::SEAMS_OK, or ::SEAMS_ERR_NULL, ::SEAMS_ERR_STATE,
 *         ::SEAMS_ERR_RANGE, ::SEAMS_ERR_BOX, ::SEAMS_ERR_ARG,
 *         ::SEAMS_ERR_NOMEM, ::SEAMS_ERR_INTERNAL.
 * @par Handle access
 * Non-const only because the call clears the handle's error slot on
 * entry, so seams_handle_error() reports this call. It reads the frame
 * and keeps no scratch on the handle.
 * @par Ownership
 * The caller owns every buffer. The library copies nothing and keeps no
 * pointer after the call returns.
 * @par Thread safety
 * Exclusive use of sys.
 * @par Example
 * @code
 * size_t n = 0;
 * seams_density_fractional(sys, 0, 4, 2, NULL, 0, 0.0, NULL, 0,
 *                          NULL, NULL, NULL, NULL, 0, &n);
 * @endcode
 */
SEAMS_API seams_status seams_density_fractional(seams_system *sys, int type, int bins, int axis,
                                                const double *direction, size_t direction_len,
                                                double origin, const double *masses,
                                                size_t masses_len, double *coordinate,
                                                double *number_density, double *mass_density,
                                                double *bin_volume, size_t cap, size_t *n_out);

/**
 * @brief Pore axis and centre, supplied or taken from the framework.
 *
 * axis 0, 1, or 2 selects that cell edge. axis -1 with direction NULL
 * selects the longest cell edge. Any other negative axis is rejected.
 * A non-NULL direction must be parallel
 * to a cell edge. centre_frac NULL places the centre at the periodic
 * centroid of the framework atoms in the perpendicular plane.
 * framework_type 0 uses every atom for that centroid.
 *
 * @param sys Frame handle.
 * @param framework_type Species used to place an automatic centre. Zero
 *        keeps every atom.
 * @param axis 0, 1, 2, or -1 for the longest edge when direction is NULL.
 *        Any other negative value is rejected.
 * @param direction Three Cartesian components parallel to a cell edge, or NULL.
 * @param direction_len 3 when direction is not NULL. Ignored otherwise.
 * @param centre_frac Three fractional coordinates, or NULL to detect the centre.
 * @param centre_len 3 when centre_frac is not NULL. Ignored otherwise.
 * @param axis_out Receives the unit Cartesian axis. Three doubles.
 * @param centre_out Receives the Cartesian centre. Three doubles.
 * @param centre_frac_out Receives the fractional centre. Three doubles.
 * @param axis_index Receives 0, 1, or 2.
 * @return ::SEAMS_OK, or ::SEAMS_ERR_NULL, ::SEAMS_ERR_STATE,
 *         ::SEAMS_ERR_RANGE, ::SEAMS_ERR_BOX, ::SEAMS_ERR_ARG,
 *         ::SEAMS_ERR_NOMEM, ::SEAMS_ERR_INTERNAL.
 * @par Handle access
 * Non-const only because the call clears the handle's error slot on
 * entry, so seams_handle_error() reports this call. It reads the frame
 * and keeps no scratch on the handle.
 * @par Ownership
 * The caller owns the three output arrays. Each holds three doubles.
 * The library keeps no pointer.
 * @par Thread safety
 * Exclusive use of sys.
 * @par Example
 * @code
 * double axis[3], centre[3], frac[3];
 * int index = -1;
 * seams_pore_frame(sys, 0, -1, NULL, 0, NULL, 0, axis, centre, frac, &index);
 * @endcode
 */
SEAMS_API seams_status seams_pore_frame(seams_system *sys, int framework_type, int axis,
                                        const double *direction, size_t direction_len,
                                        const double *centre_frac, size_t centre_len,
                                        double *axis_out, double *centre_out,
                                        double *centre_frac_out, int *axis_index);

/**
 * @brief Radial, azimuthal, and axial densities about a pore axis.
 *
 * The axis is a cell edge, as in seams_pore_frame(). r_max 0 covers the
 * perpendicular face. Radial bins run from 0 to r_max. The radial
 * coordinate is the distance to the nearest periodic image of the pore
 * axis. Within the inscribed radius of the centred face that image is
 * the face origin. Beyond the inscribed radius a corner point can be
 * closer to a neighbouring image, and the coordinate follows that
 * image. The radial shell is the annulus inside the Voronoi cell of the
 * perpendicular lattice. A circle-edge crossing that lands on a vertex
 * stays on that edge. Azimuthal bins cover [-pi, pi) in the centred
 * face, and each sector is clipped to that face. Azimuth zero points
 * along the perpendicular projection of cell edge (axis + 1) mod 3,
 * with the component parallel to the pore axis removed. The angle
 * increases towards the pore axis crossed with that direction. The axial coordinate
 * is the fractional cell coordinate, shifted by origin. A shell volume
 * is that area times the edge length. Mass density is the caller's mass
 * unit per cubic Angstrom. An amu mass gives amu per cubic Angstrom, and
 * times 1.66054 that is grams per cubic centimetre. Pass NULL columns
 * to query the three bin counts.
 *
 * @param sys Frame handle.
 * @param type Species filter. Zero keeps every atom.
 * @param radial_bins Positive radial bin count.
 * @param azimuthal_bins Positive azimuthal bin count.
 * @param axial_bins Positive axial bin count.
 * @param r_max Outer radius in Angstrom. Zero covers the face.
 * @param axis 0, 1, 2, or -1 when direction is NULL. Any other negative
 *        value is rejected.
 * @param direction Three Cartesian components parallel to a cell edge, or NULL.
 * @param direction_len 3 when direction is not NULL. Ignored otherwise.
 * @param centre_frac Three fractional coordinates, or NULL to detect the centre.
 * @param centre_len 3 when centre_frac is not NULL. Ignored otherwise.
 * @param framework_type Species for an automatic centre. Zero keeps every atom.
 * @param origin Fractional axial origin, wrapped into [0, 1).
 * @param masses One mass per stored atom, or NULL for mass 1.
 * @param masses_len Length of masses. Must equal the atom count when
 *        masses is not NULL.
 * @param radius Radial bin centres, or NULL to skip.
 * @param radial_number Radial number density, or NULL to skip.
 * @param radial_mass Radial mass density in the caller's mass units per
 *        cubic Angstrom, or NULL to skip. An amu mass gives amu per cubic
 *        Angstrom; times 1.66054 that is grams per cubic centimetre.
 * @param radial_volume Radial shell volume, or NULL to skip.
 * @param radial_cap Capacity of each non-NULL radial column.
 * @param n_radial Receives the radial bin count.
 * @param theta Azimuthal bin centres in radians, or NULL to skip.
 * @param azimuthal_number Azimuthal number density, or NULL to skip.
 * @param azimuthal_mass Azimuthal mass density in the caller's mass units
 *        per cubic Angstrom, or NULL to skip. An amu mass gives amu per
 *        cubic Angstrom; times 1.66054 that is grams per cubic centimetre.
 * @param azimuthal_volume Azimuthal sector volume, or NULL to skip.
 * @param azimuthal_cap Capacity of each non-NULL azimuthal column.
 * @param n_azimuthal Receives the azimuthal bin count.
 * @param axial Axial bin centres, or NULL to skip.
 * @param axial_number Axial number density, or NULL to skip.
 * @param axial_mass Axial mass density in the caller's mass units per
 *        cubic Angstrom, or NULL to skip. An amu mass gives amu per cubic
 *        Angstrom; times 1.66054 that is grams per cubic centimetre.
 * @param axial_volume Axial slab volume, or NULL to skip.
 * @param axial_cap Capacity of each non-NULL axial column.
 * @param n_axial Receives the axial bin count.
 * @return ::SEAMS_OK, or ::SEAMS_ERR_NULL, ::SEAMS_ERR_STATE,
 *         ::SEAMS_ERR_RANGE, ::SEAMS_ERR_BOX, ::SEAMS_ERR_ARG,
 *         ::SEAMS_ERR_NOMEM, ::SEAMS_ERR_INTERNAL.
 * @par Handle access
 * Non-const only because the call clears the handle's error slot on
 * entry, so seams_handle_error() reports this call. It reads the frame
 * and keeps no scratch on the handle.
 * @par Ownership
 * The caller owns every column. The library keeps no pointer.
 * @par Thread safety
 * Exclusive use of sys.
 * @par Example
 * @code
 * size_t nr = 0, nt = 0, ns = 0;
 * seams_pore_profiles(sys, 0, 4, 4, 4, 0.0, 2, NULL, 0, NULL, 0, 0, 0.0, NULL, 0,
 *                     NULL, NULL, NULL, NULL, 0, &nr,
 *                     NULL, NULL, NULL, NULL, 0, &nt,
 *                     NULL, NULL, NULL, NULL, 0, &ns);
 * @endcode
 */
SEAMS_API seams_status seams_pore_profiles(
    seams_system *sys, int type, int radial_bins, int azimuthal_bins, int axial_bins, double r_max,
    int axis, const double *direction, size_t direction_len, const double *centre_frac,
    size_t centre_len, int framework_type, double origin, const double *masses, size_t masses_len,
    double *radius, double *radial_number, double *radial_mass, double *radial_volume,
    size_t radial_cap, size_t *n_radial, double *theta, double *azimuthal_number,
    double *azimuthal_mass, double *azimuthal_volume, size_t azimuthal_cap, size_t *n_azimuthal,
    double *axial, double *axial_number, double *axial_mass, double *axial_volume, size_t axial_cap,
    size_t *n_axial);

/**
 * @brief Nearest framework site and contact counts for each guest.
 *
 * A lattice coefficient at or above 2^52, or a cell edge past the Euclidean
 * range of binary64, returns SEAMS_ERR_RANGE. The fractional wrap is not
 * substituted.
 *
 * guest_type must be positive. framework_type 0 keeps every atom of
 * another type. framework_type equal to guest_type is an error.
 * nearest_distance is negative when a guest has no
 * framework atom. A contact is a framework atom at or inside cutoff
 * under the triclinic minimum image. contact_counts is row-major,
 * one row per guest and one column per framework type, types sorted.
 * Pass NULL columns to query the guest and type counts.
 *
 * @param sys Frame handle.
 * @param guest_type Positive guest species.
 * @param framework_type Framework species. Zero keeps every other type.
 *        A value equal to guest_type is rejected.
 * @param cutoff Contact distance in Angstrom. Zero counts only a site
 *        at the same point. Negative is rejected.
 * @param guest_index Guest indices, or NULL to skip.
 * @param nearest_index Nearest framework index, or -1. NULL skips the column.
 * @param nearest_type Type of that site, or NULL to skip.
 * @param nearest_distance Distance in Angstrom, or NULL to skip.
 * @param guest_cap Capacity of each non-NULL guest column.
 * @param n_guest Receives the guest count.
 * @param framework_types Sorted framework types, or NULL to skip.
 * @param type_cap Capacity of framework_types when that pointer is not NULL.
 * @param n_types Receives the framework type count.
 * @param contact_counts Row-major counts, or NULL to skip.
 * @param contact_cap Capacity of contact_counts. Must hold n_guest * n_types
 *        when contact_counts is not NULL.
 * @return ::SEAMS_OK, or ::SEAMS_ERR_NULL, ::SEAMS_ERR_STATE,
 *         ::SEAMS_ERR_RANGE, ::SEAMS_ERR_BOX, ::SEAMS_ERR_ARG,
 *         ::SEAMS_ERR_NOMEM, ::SEAMS_ERR_INTERNAL.
 * @par Handle access
 * Non-const only because the call clears the handle's error slot on
 * entry, so seams_handle_error() reports this call. It reads the frame
 * and keeps no scratch on the handle.
 * @par Ownership
 * The caller owns every column. The library keeps no pointer.
 * @par Thread safety
 * Exclusive use of sys.
 * @par Example
 * @code
 * size_t ng = 0, nt = 0;
 * seams_host_guest(sys, 1, 0, 3.5, NULL, NULL, NULL, NULL, 0, &ng,
 *                  NULL, 0, &nt, NULL, 0);
 * @endcode
 */
SEAMS_API seams_status seams_host_guest(seams_system *sys, int guest_type, int framework_type,
                                        double cutoff, int *guest_index, int *nearest_index,
                                        int *nearest_type, double *nearest_distance,
                                        size_t guest_cap, size_t *n_guest, int *framework_types,
                                        size_t type_cap, size_t *n_types, int *contact_counts,
                                        size_t contact_cap);

/** @} */


#ifdef __cplusplus
}
#endif

#endif
