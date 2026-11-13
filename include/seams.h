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
 * ABI policy. SEAMS_ABI stays 1 while every shipped function keeps its
 * signature and every shipped enumerator keeps its numeric value. New
 * functions and new enumerators are appended. No public struct carries
 * fields. A later struct, if one is added, starts with a uint32_t size
 * set by the caller to the size of the struct the caller compiled, and
 * the library reads only the prefix that size covers. Status codes are
 * the seams_status enumeration. A failing call stores a detail string
 * for seams_last_error() and, when the call had a handle, for
 * seams_handle_error(). Exceptions do not cross this boundary.
 *
 * Buffer convention for functions added after the first ABI 1 set.
 * Lengths are size_t. A NULL output buffer with a non-NULL length
 * pointer queries the required count and returns SEAMS_OK. A short
 * non-NULL buffer returns SEAMS_ERR_RANGE and writes nothing into that
 * buffer. The length pointer still receives the required count.
 *
 * The version macros live in seams_version.h, which this header
 * includes: SEAMS_VERSION_MAJOR, SEAMS_VERSION_MINOR,
 * SEAMS_VERSION_PATCH, SEAMS_VERSION_STRING, SEAMS_VERSION_NUMBER.
 *
 * SPDX-License-Identifier: MIT
 * Copyright (c) d-SEAMS developers
 */

#define SEAMS_ABI 1
#include "seams_version.h"

#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

#if defined(_WIN32)
#define SEAMS_API __declspec(dllexport)
#elif defined(__GNUC__) || defined(__clang__)
#define SEAMS_API __attribute__((visibility("default")))
#else
#define SEAMS_API
#endif

#if defined(__GNUC__) || defined(__clang__)
#define SEAMS_DEPRECATED(msg) __attribute__((deprecated(msg)))
#elif defined(_MSC_VER)
#define SEAMS_DEPRECATED(msg) __declspec(deprecated(msg))
#else
#define SEAMS_DEPRECATED(msg)
#endif

/** Status codes. Zero is success. */
typedef enum seams_status {
  SEAMS_OK = 0,            /**< The call succeeded. */
  SEAMS_ERR_NULL = 1,      /**< A required pointer is null. */
  SEAMS_ERR_RANGE = 2,     /**< A number is outside the accepted range. */
  SEAMS_ERR_BOX = 3,       /**< The cell is missing or has no volume. */
  SEAMS_ERR_ARG = 4,       /**< A caller argument is the wrong shape. */
  SEAMS_ERR_STATE = 5,     /**< The handle has no frame yet. */
  SEAMS_ERR_NOMEM = 6,     /**< An allocation failed. */
  SEAMS_ERR_INTERNAL = 7   /**< The library hit an internal failure. */
} seams_status;

/** CHILL and CHILL+ atom labels. Matches molSys::atom_state_type. */
typedef enum seams_label {
  SEAMS_LABEL_CUBIC = 0,
  SEAMS_LABEL_HEXAGONAL = 1,
  SEAMS_LABEL_WATER = 2,
  SEAMS_LABEL_INTERFACIAL = 3,
  SEAMS_LABEL_CLATHRATE = 4,
  SEAMS_LABEL_INTERCLATHRATE = 5,
  SEAMS_LABEL_UNCLASSIFIED = 6,
  SEAMS_LABEL_RECUBIC = 7,
  SEAMS_LABEL_REHEX = 8
} seams_label;

/** Ring topology. Size classes are used when a six-ring is not a cage face. */
typedef enum seams_ring_symmetry {
  SEAMS_RING_UNCLASSIFIED = 0,
  SEAMS_RING_DDC = 1,
  SEAMS_RING_HC_BASAL = 2,
  SEAMS_RING_HC_PRISMATIC = 3,
  SEAMS_RING_BOTH_BASAL = 4,
  SEAMS_RING_BOTH_PRISMATIC = 5,
  SEAMS_RING_PRISM = 6,
  SEAMS_RING_DEFORMED_PRISM = 7,
  SEAMS_RING_MIXED_PRISM = 8,
  SEAMS_RING_PENTAGON = 9,
  SEAMS_RING_HEXAGON = 10,
  SEAMS_RING_OTHER = 11
} seams_ring_symmetry;

/** Published edges are inclusive. OVITO edges are open. */
typedef enum seams_chill_rule {
  SEAMS_CHILL_PUBLISHED = 0,
  SEAMS_CHILL_OVITO = 1
} seams_chill_rule;

/** Neighbour graph. Seeded is a cage graph: mutual seeds, union completions. */
typedef enum seams_graph {
  SEAMS_GRAPH_CUTOFF = 0,
  SEAMS_GRAPH_MUTUAL = 1,
  SEAMS_GRAPH_SEEDED = 2
} seams_graph;

/**
 * Per-atom cage membership from the hexagonal and double-diamond search.
 * The four codes are disjoint.
 */
typedef enum seams_cage_label {
  SEAMS_CAGE_WATER = 0, /**< In neither cage. */
  SEAMS_CAGE_HEX = 1,   /**< Hexagonal cage only. */
  SEAMS_CAGE_CUBIC = 2, /**< Double-diamond cage only. */
  SEAMS_CAGE_BOTH = 3   /**< Both cages. */
} seams_cage_label;

/** One frame. The caller creates and destroys it. */
typedef struct seams_system seams_system;
/** One neighbour list owned by the caller, built from a frame. */
typedef struct seams_neighbours seams_neighbours;
/** One primitive-ring set owned by the caller. */
typedef struct seams_rings seams_rings;
/** Hexagonal and double-diamond membership for one frame. */
typedef struct seams_cage_set seams_cage_set;
/** One planar ice/brine front measurement. */
typedef struct seams_brine seams_brine;
/** Topological fingerprint of one bonded frame. */
typedef struct seams_topo seams_topo;

/** Site kind. Values match the engine's site table. */
typedef enum seams_site_kind {
  SEAMS_SITE_UNSPECIFIED = 0,
  SEAMS_SITE_CATION = 1,
  SEAMS_SITE_ANION = 2,
  SEAMS_SITE_TAIL = 3,
  SEAMS_SITE_DONOR_H = 4,
  SEAMS_SITE_ACCEPTOR = 5,
  SEAMS_SITE_POLAR = 6,
  SEAMS_SITE_APOLAR = 7,
  SEAMS_SITE_WATER_O = 8,
  SEAMS_SITE_WATER_H = 9,
  SEAMS_SITE_SOLVENT = 10
} seams_site_kind;

/**
 * @brief ABI number compiled into this library.
 * @return SEAMS_ABI, which is 1.
 * Thread safety: pure.
 */
SEAMS_API int seams_abi(void);

/**
 * @brief Packed release number.
 * @return SEAMS_VERSION_NUMBER, (major << 16) | (minor << 8) | patch.
 * Thread safety: pure.
 */
SEAMS_API int seams_version(void);

/** @brief Major version. Thread safety: pure. */
SEAMS_API int seams_version_major(void);
/** @brief Minor version. Thread safety: pure. */
SEAMS_API int seams_version_minor(void);
/** @brief Patch version. Thread safety: pure. */
SEAMS_API int seams_version_patch(void);

/**
 * @brief Version text, for example "2.14.0".
 * @return Pointer to a static string. The caller does not free it.
 * Thread safety: pure.
 */
SEAMS_API const char *seams_version_string(void);

/**
 * @brief Short name of a status code.
 * @param status A seams_status value.
 * @return Pointer to a static string. Unknown codes return "unknown status".
 * Thread safety: pure.
 */
SEAMS_API const char *seams_status_string(seams_status status);

/**
 * @brief Thread-local detail of the last failing call on this thread.
 *
 * The string is empty after a successful call on this thread. It is
 * owned by the library and remains valid until the next call on this
 * thread that sets or clears the error. Prefer seams_handle_error()
 * when the failing call had a handle: that copy survives a later call
 * on another handle.
 *
 * @return Pointer to a thread-local string. Never NULL. The caller does
 *         not free it.
 * Thread safety: the storage is thread-local.
 */
SEAMS_API const char *seams_last_error(void);

/**
 * @brief Detail string stored on a handle.
 * @param handle A seams_system, seams_neighbours, seams_rings,
 *        seams_cage_set, seams_topo, or seams_brine, or NULL.
 * @return The handle's string, or the thread-local string when handle
 *         is NULL. Never NULL. The caller does not free it. The handle
 *         string stays valid until the next call on that handle or
 *         until destroy.
 * Thread safety: reading a handle requires that no other thread is
 * inside a call on that same handle.
 */
SEAMS_API const char *seams_handle_error(const void *handle);

/**
 * @brief Name of a CHILL label.
 * @param label A seams_label value.
 * @return Static string. Unknown codes return "unknown".
 * Thread safety: pure.
 */
SEAMS_API const char *seams_label_name(int label);

/**
 * @brief Name of a ring-symmetry code.
 * @param symmetry A seams_ring_symmetry value.
 * @return Static string. Unknown codes return "unknown".
 * Thread safety: pure.
 */
SEAMS_API const char *seams_ring_symmetry_name(int symmetry);

/**
 * @brief Set the OpenMP team used by later analysis calls.
 * @param n Positive width, or 0 to restore the width captured at first
 *        use. A negative width is rejected.
 * @return SEAMS_OK, or SEAMS_ERR_RANGE when n is negative.
 * Thread safety: process-global. Do not call it concurrently with
 * analysis. The stored width is then read by every handle.
 */
SEAMS_API seams_status seams_set_num_threads(int n);

/**
 * @brief Team width the next analysis call will request.
 * @return A positive count. After seams_set_num_threads(0) this is the
 *         saved default, which is omp_get_max_threads() when OpenMP is
 *         linked and 1 otherwise.
 * Thread safety: observes the process-global width. Do not call it
 * concurrently with seams_set_num_threads().
 */
SEAMS_API int seams_get_num_threads(void);

/**
 * @brief Allocate an empty frame.
 * @return A new handle, or NULL on allocation failure. seams_last_error()
 *         then explains the failure. The caller owns the handle and
 *         passes it to seams_system_destroy().
 * Thread safety: safe. The new handle is not shared yet.
 */
SEAMS_API seams_system *seams_system_create(void);

/**
 * @brief Free a frame.
 * @param sys Handle from seams_system_create(), or NULL.
 * Ownership: the handle and the frame it stored are released. Neighbour,
 * ring, cage, and brine handles built from it stay valid: they hold
 * their own copies.
 * Thread safety: the caller must not use sys after this returns, and
 * must not call it concurrently with any use of sys.
 */
SEAMS_API void seams_system_destroy(seams_system *sys);

/**
 * @brief Store points and a general cell. Lengths are size_t.
 * @param sys Frame handle.
 * @param xyz n times 3 Cartesian coordinates, Angstrom.
 * @param n Atom count. Zero is rejected. A count that does not fit in
 *        the frame's integer width is rejected.
 * @param box9 Cell numbers. box_len selects how many are read.
 * @param box_len Must be 9.
 * @param types Atom types, or NULL to store type 1 for every atom.
 * @param ids Atom ids, or NULL to store 1..n. Ids must be unique.
 * @return SEAMS_OK, or SEAMS_ERR_NULL, SEAMS_ERR_RANGE, SEAMS_ERR_BOX,
 *         SEAMS_ERR_ARG, SEAMS_ERR_NOMEM, SEAMS_ERR_INTERNAL.
 * Ownership: the library copies the arrays. The caller keeps them.
 * Thread safety: exclusive use of sys.
 */
SEAMS_API seams_status seams_system_set_atoms(seams_system *sys, const double *xyz,
                                              size_t n, const double *box9,
                                              size_t box_len, const int *types,
                                              const int *ids);

/**
 * @brief Store points and a general LAMMPS cell.
 *
 * ABI 1 entry. n is an int. New callers use seams_system_set_atoms().
 * @param sys Frame handle.
 * @param xyz n times 3 coordinates, Angstrom.
 * @param n Atom count. A non-positive count is rejected.
 * @param box9 Nine cell numbers: ox, oy, oz, lx, ly, lz, xy, xz, yz.
 * @param types Types, or NULL for type 1.
 * @param ids Ids, or NULL for 1..n.
 * @return Same codes as seams_system_set_atoms().
 * Ownership: the library copies the arrays.
 * Thread safety: exclusive use of sys.
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
 * Thread safety: exclusive use of sys.
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
 * @return SEAMS_OK, or SEAMS_ERR_NULL, SEAMS_ERR_RANGE, SEAMS_ERR_ARG,
 *         SEAMS_ERR_NOMEM, SEAMS_ERR_INTERNAL.
 * Ownership: the library copies the frame out of the file.
 * Thread safety: exclusive use of sys. Concurrent reads of distinct
 * paths on distinct handles are safe.
 */
SEAMS_API seams_status seams_system_read(seams_system *sys, const char *path,
                                         int frame, int type_filter);

/**
 * @brief Number of frames in a trajectory path.
 * @param path Filesystem path.
 * @param nframes Receives the count. Must be non-NULL.
 * @return SEAMS_OK, or SEAMS_ERR_NULL, SEAMS_ERR_RANGE, SEAMS_ERR_ARG,
 *         SEAMS_ERR_INTERNAL. *nframes is unchanged on failure.
 * Thread safety: safe for distinct paths.
 */
SEAMS_API seams_status seams_path_nframes(const char *path, size_t *nframes);

/**
 * @brief Number of points in the frame.
 * @param sys Frame handle, or NULL.
 * @return The atom count, or 0 when sys is NULL or holds no frame.
 * Thread safety: exclusive use of sys.
 */
SEAMS_API int seams_system_count(const seams_system *sys);

/**
 * @brief Build a cutoff or mutual neighbour list.
 * @param sys Frame handle that already holds a frame.
 * @param cutoff Positive cutoff, Angstrom. The mutual graph keeps the
 *        four nearest neighbours inside cutoff + 1.5.
 * @param type Species to keep. Zero and negative keep every atom.
 * @param graph SEAMS_GRAPH_CUTOFF or SEAMS_GRAPH_MUTUAL.
 *        SEAMS_GRAPH_SEEDED is rejected here.
 * @param out Receives a new handle. On failure *out is NULL.
 * @return SEAMS_OK, or SEAMS_ERR_NULL, SEAMS_ERR_STATE, SEAMS_ERR_RANGE,
 *         SEAMS_ERR_ARG, SEAMS_ERR_NOMEM, SEAMS_ERR_INTERNAL.
 * Ownership: the caller frees *out with seams_neighbours_destroy().
 * The list is a copy and does not borrow sys.
 * Thread safety: exclusive use of sys. The new handle is private to
 * the caller.
 */
SEAMS_API seams_status seams_neighbours_build(seams_system *sys, double cutoff,
                                              int type, int graph,
                                              seams_neighbours **out);

/**
 * @brief Free a neighbour list.
 * @param nb Handle from seams_neighbours_build(), or NULL.
 * Thread safety: the caller must not use nb after this returns.
 */
SEAMS_API void seams_neighbours_destroy(seams_neighbours *nb);

/**
 * @brief Number of neighbour rows.
 * @param nb Neighbour handle, or NULL.
 * @return The row count, or 0 when nb is NULL.
 * Thread safety: exclusive use of nb.
 */
SEAMS_API int seams_neighbours_count(const seams_neighbours *nb);

/**
 * @brief Borrow the neighbour list as CSR.
 *
 * Offsets has n + 1 entries. Columns omit the self entry. The pointers
 * are owned by nb and remain valid until seams_neighbours_destroy(nb).
 * New callers that need their own buffers use seams_neighbours_copy().
 *
 * @param nb Neighbour handle.
 * @param offsets Receives the offset pointer.
 * @param n_offsets Receives the offset count, n + 1.
 * @param cols Receives the column pointer. May be NULL when the list
 *        has no pairs; *nnz is then 0.
 * @param nnz Receives the column count.
 * @return SEAMS_OK, or SEAMS_ERR_NULL.
 * Thread safety: exclusive use of nb. The returned pointers must not
 * be used concurrently with a call that frees nb.
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
 * @return SEAMS_OK when every non-NULL buffer was filled, or when both
 *         buffers were NULL and the lengths were written.
 *         SEAMS_ERR_RANGE when a non-NULL buffer is short.
 *         SEAMS_ERR_NULL when nb, n_offsets, or nnz is NULL.
 * Ownership: the caller owns the buffers. The handle is not modified.
 * Thread safety: exclusive use of nb.
 */
SEAMS_API seams_status seams_neighbours_copy(const seams_neighbours *nb,
                                             int *offsets, size_t offsets_cap,
                                             size_t *n_offsets, int *cols,
                                             size_t cols_cap, size_t *nnz);

/**
 * @brief Original CHILL labels, published or OVITO edges.
 *
 * Published CHILL is Moore's four-neighbour rule with inclusive c_ij
 * edges. OVITO mode uses the same thresholds and the open edges.
 * Labels match seams_label. counts9, when non-NULL, receives nine
 * population counts in label order.
 *
 * @param sys Frame handle.
 * @param cutoff Positive cutoff, Angstrom.
 * @param type Species filter. Zero and negative keep every atom.
 * @param graph SEAMS_GRAPH_CUTOFF or SEAMS_GRAPH_MUTUAL.
 * @param rule SEAMS_CHILL_PUBLISHED or SEAMS_CHILL_OVITO.
 * @param labels Caller buffer of labels_len ints.
 * @param labels_len Must be at least the frame count.
 * @param counts9 Nine ints, or NULL to skip the histogram.
 * @return SEAMS_OK, or SEAMS_ERR_NULL, SEAMS_ERR_STATE, SEAMS_ERR_RANGE,
 *         SEAMS_ERR_ARG, SEAMS_ERR_NOMEM, SEAMS_ERR_INTERNAL.
 * Ownership: the caller owns labels and counts9.
 * Thread safety: exclusive use of sys.
 */
SEAMS_API seams_status seams_chill(seams_system *sys, double cutoff, int type,
                                   int graph, int rule, int *labels,
                                   size_t labels_len, int counts9[9]);

/**
 * @brief CHILL+ labels.
 *
 * rule selects the published Nguyen and Molinero table or the OVITO
 * Chill+ modifier. counts9, when non-NULL, receives nine population
 * counts. labels has one entry per atom in the frame; the length is
 * the frame count.
 *
 * @param sys Frame handle.
 * @param cutoff Positive cutoff, Angstrom.
 * @param type Species filter. Zero and negative keep every atom.
 * @param graph SEAMS_GRAPH_CUTOFF or SEAMS_GRAPH_MUTUAL.
 * @param rule SEAMS_CHILL_PUBLISHED or SEAMS_CHILL_OVITO.
 * @param labels Caller buffer of one label per atom.
 * @param counts9 Nine ints, or NULL.
 * @return SEAMS_OK, or SEAMS_ERR_NULL, SEAMS_ERR_STATE, SEAMS_ERR_RANGE,
 *         SEAMS_ERR_ARG, SEAMS_ERR_NOMEM, SEAMS_ERR_INTERNAL.
 * Ownership: the caller owns labels and counts9.
 * Thread safety: exclusive use of sys.
 */
SEAMS_API seams_status seams_chill_plus(seams_system *sys, double cutoff, int type,
                                        int graph, int rule, int *labels,
                                        int counts9[9]);

/**
 * @brief Local Steinhardt q_l, and the Lechner-Dellago average on request.
 *
 * order_l is 3, 4, 6, 8, or 12. Orders 3, 4, and 6 are q3, q4, and q6.
 * Either array may be NULL. A non-NULL array has one double per atom
 * and its length must be at least the frame count. The average is the
 * neighbour mean of q_lm. The cutoff graph reads the cutoff shell. The
 * mutual graph reads the four nearest neighbours inside cutoff + 1.5.
 *
 * @param sys Frame handle.
 * @param cutoff Positive cutoff, Angstrom.
 * @param type Species filter. Zero and negative keep every atom.
 * @param order_l Spherical-harmonic degree.
 * @param graph SEAMS_GRAPH_CUTOFF or SEAMS_GRAPH_MUTUAL.
 * @param ql Local q_l, or NULL.
 * @param ql_len Length of ql. Ignored when ql is NULL.
 * @param qlbar Averaged q_l, or NULL to skip the average.
 * @param qlbar_len Length of qlbar. Ignored when qlbar is NULL.
 * @return SEAMS_OK, or SEAMS_ERR_NULL, SEAMS_ERR_STATE, SEAMS_ERR_RANGE,
 *         SEAMS_ERR_ARG, SEAMS_ERR_NOMEM, SEAMS_ERR_INTERNAL.
 * Ownership: the caller owns both arrays.
 * Thread safety: exclusive use of sys. The spherical-harmonic scratch
 * is process-global, so this call does not run beside seams_q4q6() or
 * seams_q6bar().
 */
SEAMS_API seams_status seams_ql(seams_system *sys, double cutoff, int type,
                                int order_l, int graph, double *ql, size_t ql_len,
                                double *qlbar, size_t qlbar_len);

/**
 * @brief Local q4 and q6 on the cutoff shell.
 *
 * ABI 1 pair. Either pointer may be NULL, and the two arrays may not
 * both be NULL. Each non-NULL array has one double per atom. Degrees
 * 4 and 6 use the solid-harmonic kernel. The Lechner-Dellago averages
 * are seams_ql() with a non-NULL qlbar, and the nucleation summary is
 * seams_q6bar().
 *
 * @param sys Frame handle.
 * @param cutoff Positive cutoff, Angstrom.
 * @param type Species filter.
 * @param q4 Local q4, or NULL.
 * @param q6 Local q6, or NULL.
 * @return SEAMS_OK, or SEAMS_ERR_NULL, SEAMS_ERR_STATE, SEAMS_ERR_RANGE,
 *         SEAMS_ERR_NOMEM, SEAMS_ERR_INTERNAL.
 * Ownership: the caller owns q4 and q6.
 * Thread safety: exclusive use of sys. Shares the process-global
 * scratch with seams_ql() and seams_q6bar().
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
 * @param cutoff Positive cutoff, Angstrom.
 * @param type Species filter.
 * @param ice_cut Finite threshold on q6-bar.
 * @param q6bar Per-atom average, or NULL to keep only the summary.
 * @param q6bar_len Length of q6bar. Ignored when q6bar is NULL.
 * @param mean Receives the mean, or NULL.
 * @param median Receives the median, or NULL.
 * @param n_ice Receives the count above ice_cut, or NULL.
 * @return SEAMS_OK, or SEAMS_ERR_NULL, SEAMS_ERR_STATE, SEAMS_ERR_RANGE,
 *         SEAMS_ERR_NOMEM, SEAMS_ERR_INTERNAL. At least one of q6bar,
 *         mean, median, and n_ice must be non-NULL.
 * Ownership: the caller owns q6bar.
 * Thread safety: exclusive use of sys. Shares the process-global
 * scratch with seams_ql() and seams_q4q6().
 */
SEAMS_API seams_status seams_q6bar(seams_system *sys, double cutoff, int type,
                                   double ice_cut, double *q6bar, size_t q6bar_len,
                                   double *mean, double *median, int *n_ice);

/**
 * @brief Primitive rings (Franzblau) up to a maximum size.
 * @param sys Frame handle.
 * @param cutoff Positive cutoff, Angstrom.
 * @param type Species filter.
 * @param graph SEAMS_GRAPH_CUTOFF or SEAMS_GRAPH_MUTUAL.
 * @param max_depth Largest ring size, from 3 through 12.
 * @param out Receives a new handle. On failure *out is NULL.
 * @return SEAMS_OK, or SEAMS_ERR_NULL, SEAMS_ERR_STATE, SEAMS_ERR_RANGE,
 *         SEAMS_ERR_ARG, SEAMS_ERR_NOMEM, SEAMS_ERR_INTERNAL.
 * Ownership: the caller frees *out with seams_rings_destroy().
 * Thread safety: exclusive use of sys.
 */
SEAMS_API seams_status seams_rings_find(seams_system *sys, double cutoff, int type,
                                        int graph, int max_depth, seams_rings **out);

/**
 * @brief Free a ring set.
 * @param rings Handle from seams_rings_find(), or NULL.
 * Thread safety: the caller must not use rings after this returns.
 */
SEAMS_API void seams_rings_destroy(seams_rings *rings);

/**
 * @brief Number of primitive rings.
 * @param rings Ring handle, or NULL.
 * @return The ring count, or 0 when rings is NULL.
 * Thread safety: exclusive use of rings.
 */
SEAMS_API int seams_rings_count(const seams_rings *rings);

/**
 * @brief Count primitive rings by size.
 * @param rings Ring handle.
 * @param counts Caller buffer. counts[k] receives the number of rings
 *        of size k. Sizes outside the buffer are ignored.
 * @param counts_len Length of counts. Must be positive.
 * @return SEAMS_OK, or SEAMS_ERR_NULL, SEAMS_ERR_RANGE.
 * Ownership: the caller owns counts.
 * Thread safety: exclusive use of rings.
 */
SEAMS_API seams_status seams_rings_histogram(const seams_rings *rings, int *counts,
                                             int counts_len);

/**
 * @brief Size and symmetry of each ring.
 * @param rings Ring handle.
 * @param sizes One int per ring, or NULL.
 * @param symmetry One seams_ring_symmetry code per ring, or NULL.
 *        The two pointers may not both be NULL. Each non-NULL buffer
 *        has length seams_rings_count(rings).
 * @return SEAMS_OK, or SEAMS_ERR_NULL.
 * Ownership: the caller owns both arrays.
 * Thread safety: exclusive use of rings.
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
 * @return SEAMS_OK, or SEAMS_ERR_NULL, SEAMS_ERR_RANGE.
 * Ownership: the caller owns atoms. Indices are positions in the frame.
 * Thread safety: exclusive use of rings.
 */
SEAMS_API seams_status seams_rings_vertices(const seams_rings *rings, int index,
                                            int *atoms, size_t cap, size_t *n_out);

/**
 * @brief Fraction of atoms that sit on a ring.
 * @param rings Ring handle.
 * @param size Ring size to count. Zero counts atoms that sit on any ring.
 * @param fraction Receives a value in [0, 1].
 * @return SEAMS_OK, or SEAMS_ERR_NULL.
 * Thread safety: exclusive use of rings.
 */
SEAMS_API seams_status seams_rings_fraction(const seams_rings *rings, int size,
                                            double *fraction);

/**
 * @brief Hexagonal and double-diamond cage membership.
 * @param sys Frame handle.
 * @param cutoff Positive cutoff, Angstrom. The seeded graph uses
 *        cutoff + 1.5 as the candidate shell.
 * @param type Species filter.
 * @param graph SEAMS_GRAPH_CUTOFF, SEAMS_GRAPH_MUTUAL, or
 *        SEAMS_GRAPH_SEEDED.
 * @param k Neighbour count for the seeded graph. Zero selects 4.
 * @param complete Non-zero runs ring-adjacent completion on the seeded
 *        labels. Ignored for the cutoff and mutual graphs.
 * @param out Receives a new handle. On failure *out is NULL.
 * @return SEAMS_OK, or SEAMS_ERR_NULL, SEAMS_ERR_STATE, SEAMS_ERR_RANGE,
 *         SEAMS_ERR_ARG, SEAMS_ERR_NOMEM, SEAMS_ERR_INTERNAL.
 * Ownership: the caller frees *out with seams_cages_destroy().
 * Thread safety: exclusive use of sys.
 */
SEAMS_API seams_status seams_cages_find(seams_system *sys, double cutoff, int type,
                                        int graph, int k, int complete,
                                        seams_cage_set **out);

/**
 * @brief Free a cage-membership handle.
 * @param cages Handle from seams_cages_find(), or NULL.
 * Thread safety: the caller must not use cages after this returns.
 */
SEAMS_API void seams_cages_destroy(seams_cage_set *cages);

/**
 * @brief Disjoint cage population counts.
 * @param cages Cage handle.
 * @param n_hex Hexagonal only, or NULL.
 * @param n_cubic Double-diamond only, or NULL.
 * @param n_water Neither cage, or NULL.
 * @param n_both Both cages, or NULL.
 * @return SEAMS_OK, or SEAMS_ERR_NULL when cages is NULL. At least one
 *         count pointer must be non-NULL.
 * The four counts sum to the number of atoms the search kept.
 * Thread safety: exclusive use of cages.
 */
SEAMS_API seams_status seams_cages_counts(const seams_cage_set *cages, int *n_hex,
                                          int *n_cubic, int *n_water, int *n_both);

/**
 * @brief Per-atom cage code.
 * @param cages Cage handle.
 * @param per_atom Caller buffer of seams_cage_label values, one per atom
 *        in the frame. Atoms the species filter dropped are
 *        SEAMS_CAGE_WATER.
 * @param len Length of per_atom. Must be at least the frame count stored
 *        on the handle.
 * @return SEAMS_OK, or SEAMS_ERR_NULL, SEAMS_ERR_RANGE.
 * Ownership: the caller owns per_atom.
 * Thread safety: exclusive use of cages.
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
 * @param cutoff Positive cutoff, Angstrom.
 * @param type Species filter.
 * @param graph SEAMS_GRAPH_CUTOFF or SEAMS_GRAPH_MUTUAL.
 * @param spec Signature text. Must be non-NULL.
 * @param n_cages Receives the cage count.
 * @param n_atoms Receives the vertex count.
 * @return SEAMS_OK, or SEAMS_ERR_NULL, SEAMS_ERR_STATE, SEAMS_ERR_RANGE,
 *         SEAMS_ERR_ARG, SEAMS_ERR_NOMEM, SEAMS_ERR_INTERNAL.
 * Thread safety: exclusive use of sys.
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
 * @return SEAMS_OK, or SEAMS_ERR_NULL, SEAMS_ERR_STATE, SEAMS_ERR_RANGE,
 *         SEAMS_ERR_INTERNAL, SEAMS_ERR_NOMEM.
 * Ownership: the caller owns labels.
 * Thread safety: exclusive use of sys.
 */
SEAMS_API seams_status seams_system_symmetry(seams_system *sys, double shell_rmsd,
                                             int *labels, size_t labels_len);

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
 * @return SEAMS_OK, or SEAMS_ERR_NULL, SEAMS_ERR_STATE, SEAMS_ERR_RANGE,
 *         SEAMS_ERR_ARG, SEAMS_ERR_NOMEM, SEAMS_ERR_INTERNAL.
 * Ownership: the caller frees *out with seams_brine_destroy().
 * Thread safety: exclusive use of sys.
 */
SEAMS_API seams_status seams_brine_analyse(seams_system *sys, const int *ice,
                                           size_t n_ice, int water_type, int ion_type,
                                           int n_bins, int axis, seams_brine **out);

/**
 * @brief Free a brine-front handle.
 * @param front Handle from seams_brine_analyse(), or NULL.
 * Thread safety: the caller must not use front after this returns.
 */
SEAMS_API void seams_brine_destroy(seams_brine *front);

/**
 * @brief The two planar crossings of the ice slab.
 * @param front Brine handle.
 * @param z_enter Where the ice fraction crosses the cut on entering
 *        the ice along the positive axis, or NULL.
 * @param z_exit The crossing on leaving the ice, or NULL.
 * @param thickness Ice run along the axis, or NULL.
 * @param found Receives 1 when both crossings exist, otherwise 0, or NULL.
 * @return SEAMS_OK, or SEAMS_ERR_NULL when front is NULL. At least one
 *         output pointer must be non-NULL.
 * Lengths are Angstrom, in the frame's cell.
 * Thread safety: exclusive use of front.
 */
SEAMS_API seams_status seams_brine_fronts(const seams_brine *front, double *z_enter,
                                          double *z_exit, double *thickness, int *found);

/**
 * @brief Molal rejection coefficient k = c_ice / c_liquid.
 * @param front Brine handle.
 * @param c_ice Bulk ice molality, or NULL.
 * @param c_liquid Bulk liquid molality, or NULL.
 * @param k The ratio, or NULL. Unchanged when the ratio is undefined.
 * @param defined Receives 1 when the liquid bulk contains salt or the
 *        ice bulk does not, otherwise 0, or NULL.
 * @return SEAMS_OK, or SEAMS_ERR_NULL when front is NULL. At least one
 *         output pointer must be non-NULL.
 * Molality is formula units per kilogram of water.
 * Thread safety: exclusive use of front.
 */
SEAMS_API seams_status seams_brine_rejection(const seams_brine *front, double *c_ice,
                                             double *c_liquid, double *k, int *defined);

/**
 * @brief Concentration profile along the front axis.
 * @param front Brine handle.
 * @param z_mid Bin centres, or NULL to skip the column.
 * @param molality Bin molalities, or NULL to skip the column.
 * @param ice_fraction Bin ice fractions, or NULL to skip the column.
 * @param cap Capacity of each non-NULL column.
 * @param n_out Receives the bin count.
 * @return SEAMS_OK when every non-NULL column was filled, or when every
 *         column pointer was NULL and n_out was written.
 *         SEAMS_ERR_RANGE when a non-NULL column is short.
 *         SEAMS_ERR_NULL when front or n_out is NULL.
 * z_mid is Angstrom. molality is mol/kg. ice_fraction is dimensionless.
 * Ownership: the caller owns the columns.
 * Thread safety: exclusive use of front.
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
 * @return SEAMS_OK, or SEAMS_ERR_NULL when front is NULL. At least one
 *         output pointer must be non-NULL.
 * Thread safety: exclusive use of front.
 */
SEAMS_API seams_status seams_brine_pockets(const seams_brine *front, int *n_pockets,
                                           int *reservoir_waters);

/**
 * @brief Hexagonal and double-diamond cages. Same search as seams_cages_find().
 * @param out Receives a seams_cage_set. The caller frees it with
 *        seams_cages_destroy().
 * Thread safety: exclusive use of sys.
 */
SEAMS_API seams_status seams_cages(seams_system *sys, double cutoff, int type, int graph,
                                   int k, int complete, seams_cage_set **out);

/**
 * @brief Closed cages of one ring census. Same search as seams_cages_signature().
 * Thread safety: exclusive use of sys.
 */
SEAMS_API seams_status seams_signature_cages(seams_system *sys, double cutoff, int type,
                                             int graph, const char *spec, size_t *n_cages,
                                             size_t *n_atoms);

/**
 * @brief Local Steinhardt q_l on the cutoff shell, and the average on request.
 *
 * order_l is 3, 4, 6, 8, or 12. qlbar NULL skips the Lechner-Dellago average.
 * This is seams_ql() with the cutoff graph.
 * Thread safety: exclusive use of sys. Shares the process-global scratch
 * with seams_ql(), seams_q4q6(), and seams_q6bar().
 */
SEAMS_API seams_status seams_steinhardt(seams_system *sys, double cutoff, int type,
                                        int order_l, double *ql, size_t ql_len,
                                        double *qlbar, size_t qlbar_len);

/**
 * @brief k-nearest neighbour list.
 * @param k Positive neighbour count.
 * @param cutoff Nominal cutoff, Angstrom. Candidates are taken inside
 *        cutoff + 1.5.
 * @param mutual Non-zero keeps the mutual graph. Zero keeps the one-way
 *        k-nearest graph.
 * @param out Receives a neighbour handle. Free it with
 *        seams_neighbours_destroy().
 * Thread safety: exclusive use of sys.
 */
SEAMS_API seams_status seams_neighbours_knn(seams_system *sys, int k, double cutoff,
                                            int type, int mutual, seams_neighbours **out);

/**
 * @brief Partial radial distribution.
 * @param bins Bin count. Zero selects about rmax / 0.1.
 * @param r, g, count Caller columns, or NULL to query.
 * @param n_out Receives the bin count.
 * @param volume Receives the cell volume, or NULL.
 * r is Angstrom. g is dimensionless. volume is cubic Angstrom.
 * Thread safety: exclusive use of sys.
 */
SEAMS_API seams_status seams_rdf(seams_system *sys, double rmax, int bins, int type_i,
                                 int type_j, double *r, double *g, int *count, size_t cap,
                                 size_t *n_out, double *volume);

/**
 * @brief Running coordination number along the partial RDF.
 * Thread safety: exclusive use of sys.
 */
SEAMS_API seams_status seams_running_cn(seams_system *sys, double rmax, int bins, int type_i,
                                        int type_j, double *r, double *cn, size_t cap,
                                        size_t *n_out);

/**
 * @brief Coordination number integrated to rmax.
 * Thread safety: exclusive use of sys.
 */
SEAMS_API seams_status seams_coordination(seams_system *sys, double rmax, int bins,
                                          int type_i, int type_j, double *cn);

/**
 * @brief Hydrogen-bond count on one frame.
 * @param oxygen_type Heavy-atom species. Zero keeps every atom that is
 *        not hydrogen_type.
 * @param hydrogen_type Positive hydrogen species.
 * @param dist_cutoff Acceptor-H distance, Angstrom.
 * @param angle_cutoff O-O-H angle, degrees.
 * @param n_bonds Receives the number of hydrogen bonds.
 * Thread safety: exclusive use of sys.
 */
SEAMS_API seams_status seams_hbonds(seams_system *sys, double cutoff, int oxygen_type,
                                    int hydrogen_type, double dist_cutoff,
                                    double angle_cutoff, int *n_bonds);

/**
 * @brief Number density along one cell axis.
 * @param bins Bin count. Zero selects about one bin per 0.1 Angstrom.
 * @param axis 0, 1, or 2.
 * z is Angstrom. rho is number density.
 * Thread safety: exclusive use of sys.
 */
SEAMS_API seams_status seams_density_axis(seams_system *sys, int type, int bins, int axis,
                                          double *z, double *rho, size_t cap, size_t *n_out);

/**
 * @brief Topological fingerprint of the bonded graph.
 * @param hops Neighbourhood depth. Zero selects 2.
 * @param colour_types Non-zero colours vertices by LAMMPS type.
 * @param out Receives a seams_topo. Free it with seams_fingerprint_destroy().
 * Thread safety: exclusive use of sys.
 */
SEAMS_API seams_status seams_fingerprint(seams_system *sys, double cutoff, int type, int k,
                                         int graph, int hops, int colour_types,
                                         seams_topo **out);

/** @brief Free a fingerprint. NULL is ignored. */
SEAMS_API void seams_fingerprint_destroy(seams_topo *fp);

/**
 * @brief Frame key text. A null buffer queries *n_out, the length without
 *        the terminating NUL. A short buffer returns SEAMS_ERR_RANGE.
 */
SEAMS_API seams_status seams_fingerprint_key(const seams_topo *fp, char *buf, size_t cap,
                                             size_t *n_out);

/** @brief Local key of one atom. Same buffer rule as seams_fingerprint_key(). */
SEAMS_API seams_status seams_fingerprint_atom(const seams_topo *fp, int index, char *buf,
                                              size_t cap, size_t *n_out);

/** @brief Primitive-ring census stored on the fingerprint. counts[k] is size k. */
SEAMS_API seams_status seams_fingerprint_rings(const seams_topo *fp, int *counts,
                                               size_t counts_len);

/**
 * @brief Mutual nearest unlike contacts of a site table.
 * @param site_spec Text such as "1=cationHead,2=anion".
 * Cation and anion counts are molecules in the ion cloud.
 * Thread safety: exclusive use of sys.
 */
SEAMS_API seams_status seams_contact_pairs(seams_system *sys, const char *site_spec,
                                           size_t *n_pairs, int *n_cation, int *n_anion);

/**
 * @brief Largest domain of a polar or apolar site subset.
 * @param subset 0 is apolar, 1 is polar.
 * Thread safety: exclusive use of sys.
 */
SEAMS_API seams_status seams_domains(seams_system *sys, double cutoff, const char *site_spec,
                                     int subset, int *n_subset, int *largest,
                                     double *percolation);

/**
 * @brief Perfect and imperfect prism counts by ring size.
 * @param perfect, imperfect Length at least 7. Index is the ring size.
 *        Sizes 3 through 6 are filled.
 * Thread safety: exclusive use of sys.
 */
SEAMS_API seams_status seams_prisms(seams_system *sys, double cutoff, int type, int *perfect,
                                    int *imperfect, size_t len);

/**
 * @brief Primitive-ring histogram of a monolayer.
 * @param sheet_area Sheet area, square Angstrom. Must be positive.
 * @param counts counts[k] is the number of rings of size k.
 *        counts_len must be at least max_depth + 1.
 * Thread safety: exclusive use of sys.
 */
SEAMS_API seams_status seams_monolayer(seams_system *sys, double cutoff, int type,
                                       int max_depth, double sheet_area, int *counts,
                                       size_t counts_len);

/**
 * @brief Per-atom site kind from a site-table spec.
 * @param kinds One seams_site_kind per atom.
 * Thread safety: exclusive use of sys.
 */
SEAMS_API seams_status seams_site_table(seams_system *sys, const char *spec, int *kinds,
                                        size_t len);

#ifdef __cplusplus
}
#endif

#endif
