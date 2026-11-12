#ifndef SEAMS_H_
#define SEAMS_H_

/*
 * d-SEAMS stable C API, ABI 1.
 *
 * SPDX-License-Identifier: MIT
 * Copyright (c) d-SEAMS developers
 *
 * Opaque handles. A frame is n Cartesian points and a LAMMPS cell:
 * ox, oy, oz, lx, ly, lz, xy, xz, yz. Lengths are the edge lengths.
 * A zero tilt is orthorhombic. Coordinates are Angstrom.
 */

#define SEAMS_ABI 1
#include "seams_version.h"

#ifdef __cplusplus
extern "C" {
#endif

#if defined(_WIN32)
#define SEAMS_API __declspec(dllexport)
#else
#define SEAMS_API __attribute__((visibility("default")))
#endif

typedef enum seams_status {
  SEAMS_OK = 0,
  SEAMS_ERR_NULL = 1,
  SEAMS_ERR_RANGE = 2,
  SEAMS_ERR_BOX = 3,
  SEAMS_ERR_ARG = 4,
  SEAMS_ERR_STATE = 5,
  SEAMS_ERR_NOMEM = 6,
  SEAMS_ERR_INTERNAL = 7
} seams_status;

/* Matches molSys::atom_state_type. */
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

/* Ring topology. Size classes are used when a six-ring is not a cage face. */
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

typedef enum seams_chill_rule {
  SEAMS_CHILL_PUBLISHED = 0,
  SEAMS_CHILL_OVITO = 1
} seams_chill_rule;

typedef enum seams_graph {
  SEAMS_GRAPH_CUTOFF = 0,
  SEAMS_GRAPH_MUTUAL = 1
} seams_graph;

typedef struct seams_system seams_system;
typedef struct seams_neighbours seams_neighbours;
typedef struct seams_rings seams_rings;

SEAMS_API int seams_abi(void);
SEAMS_API int seams_version_major(void);
SEAMS_API int seams_version_minor(void);
SEAMS_API int seams_version_patch(void);
SEAMS_API const char *seams_version_string(void);
SEAMS_API const char *seams_status_string(seams_status status);
SEAMS_API const char *seams_last_error(void);
SEAMS_API const char *seams_label_name(int label);
SEAMS_API const char *seams_ring_symmetry_name(int symmetry);

/* n == 0 restores the OpenMP default (or 1). n < 0 is an error. */
SEAMS_API seams_status seams_set_num_threads(int n);
SEAMS_API int seams_get_num_threads(void);

SEAMS_API seams_system *seams_system_create(void);
SEAMS_API void seams_system_destroy(seams_system *sys);

/* types and ids may be NULL (type 1, ids 1..n). */
SEAMS_API seams_status seams_system_set_frame(seams_system *sys, const double *xyz,
                                              int n, const double *box9,
                                              const int *types, const int *ids);

SEAMS_API seams_status seams_system_set_ortho(seams_system *sys, const double *xyz,
                                              int n, double ox, double oy, double oz,
                                              double lx, double ly, double lz,
                                              const int *types, const int *ids);

SEAMS_API int seams_system_count(const seams_system *sys);

/* type <= 0 keeps every atom. graph is seams_graph. */
SEAMS_API seams_status seams_neighbours_build(seams_system *sys, double cutoff,
                                              int type, int graph,
                                              seams_neighbours **out);
SEAMS_API void seams_neighbours_destroy(seams_neighbours *nb);
SEAMS_API int seams_neighbours_count(const seams_neighbours *nb);
/* CSR without the self entry. Pointers live until destroy. */
SEAMS_API seams_status seams_neighbours_csr(const seams_neighbours *nb,
                                            const int **offsets, int *n_offsets,
                                            const int **cols, int *nnz);

/* labels length is the frame count. counts9 may be NULL. */
SEAMS_API seams_status seams_chill_plus(seams_system *sys, double cutoff, int type,
                                        int graph, int rule, int *labels,
                                        int counts9[9]);

/* Either array may be NULL. Length is the frame count.
 * Cutoff-shell q4 and q6. Degrees 4 and 6 use the solid-harmonic kernel. */
SEAMS_API seams_status seams_q4q6(seams_system *sys, double cutoff, int type,
                                  double *q4, double *q6);

/* max_depth >= 3. Primitive rings (Franzblau). */
SEAMS_API seams_status seams_rings_find(seams_system *sys, double cutoff, int type,
                                        int graph, int max_depth, seams_rings **out);
SEAMS_API void seams_rings_destroy(seams_rings *rings);
SEAMS_API int seams_rings_count(const seams_rings *rings);
SEAMS_API seams_status seams_rings_histogram(const seams_rings *rings, int *counts,
                                             int counts_len);
/* sizes and symmetry each have length seams_rings_count. Either may be NULL. */
SEAMS_API seams_status seams_rings_labels(const seams_rings *rings, int *sizes,
                                          int *symmetry);
/* size 0 counts atoms that sit on any ring. */
SEAMS_API seams_status seams_rings_fraction(const seams_rings *rings, int size,
                                            double *fraction);

#ifdef __cplusplus
}
#endif

#endif
