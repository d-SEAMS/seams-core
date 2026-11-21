#ifndef SEAMS_LANG_H_
#define SEAMS_LANG_H_

#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @file seams_lang.h
 * @brief C ABI shared by the Python, R, and Julia fronts.
 *
 * @defgroup lang_frame Frames
 * @defgroup lang_order Order parameters
 * @defgroup lang_rings Rings
 * @defgroup lang_cages Cages
 * @defgroup lang_more Further analyses
 * @defgroup lang_errors Errors
 *
 * Ice Ih (`input/traj/ice1h.lammpstrj`, every atom) is 128 hexagonal
 * waters and 256 primitive six-rings. docs/examples/lang_ice.c asserts
 * both counts. meson test docs_lang_ice compiles and runs it.
 *
 * @example lang_ice.c
 *
 * C ABI shared by the Python, R, and Julia fronts.
 * The classifiers are the engine's. This file does not implement them.
 *
 * Ice labels, same order as molSys::atom_state_type:
 *   0 cubic, 1 hexagonal, 2 water, 3 interfacial, 4 clathrate,
 *   5 interClathrate, 6 unclassified, 7 reCubic, 8 reHex.
 * Cage codes on an atom: 0 water, 1 hexagonal only, 2 double-diamond only,
 * 3 both. seams_frame_cages writes the same four disjoint counters
 * (n_hex, n_cubic, n_water, n_both); they sum to the number of kept atoms.
 * Inclusive totals are n_hex + n_both and n_cubic + n_both.
 * A box is nine doubles: ox, oy, oz, lx, ly, lz, xy, xz, yz.
 * type_filter: 0 guesses oxygen (LAMMPS type 2, else 1); >0 keeps that
 * type; <0 keeps every atom.
 * graph is "mutual" (also "knn"), "cutoff", "knn-union", or "seeded"
 * (cages only). CHILL and CHILL+ take "cutoff" (the default, used when
 * graph is NULL or "") or "mutual", the mutual four-nearest list. Mutual
 * and union bonds are retained only inside cutoff + 1.5. Returns 0 on
 * success. */

typedef struct seams_frame seams_frame;

/** Name filters and explicit type assignments for chemfiles-backed input.
 * Each pointer names an array of count UTF-8 strings. Null arrays require
 * a zero count. type_map entries use NAME=ID with a positive integer ID.
 * Include lists are ORed within a category; exclusions take precedence. */
typedef struct seams_chemfiles_load_options {
  const char *const *resnames;
  int resname_count;
  const char *const *excluded_resnames;
  int excluded_resname_count;
  const char *const *atomnames;
  int atomname_count;
  const char *const *excluded_atomnames;
  int excluded_atomname_count;
  const char *const *type_map;
  int type_map_count;
} seams_chemfiles_load_options;

/** @brief seams_lang_version. @ingroup lang_errors */
const char *seams_lang_version(void);
/** @brief seams_last_error. @ingroup lang_errors */
const char *seams_last_error(void);

/** @brief seams_path_nframes. @ingroup lang_frame */
int seams_path_nframes(const char *path);

/** Null when path is null or the file cannot be read.
 * seams_last_error() then names the reason. @ingroup lang_frame */
seams_frame *seams_frame_read(const char *path, int frame, int type_filter);
/** Read with chemfiles name filters and type assignments. options may be
 * null, which is identical to seams_frame_read. @ingroup lang_frame */
seams_frame *seams_frame_read_options(
    const char *path, int frame, int type_filter,
    const seams_chemfiles_load_options *options);
/** @brief seams_frame. @ingroup lang_frame */
seams_frame *seams_frame_from_arrays(const double *xyz, int n,
                                    const double *box9, const int *types,
                                    const int *ids);
/** @brief seams_frame_free. @ingroup lang_frame */
void seams_frame_free(seams_frame *frame);

/** @brief seams_frame_nop. @ingroup lang_frame */
int seams_frame_nop(const seams_frame *frame);
/** @brief seams_frame_index. @ingroup lang_frame */
int seams_frame_index(const seams_frame *frame);
/** @brief seams_frame_box. @ingroup lang_frame */
void seams_frame_box(const seams_frame *frame, double box9[9]);
/** @brief seams_frame_positions. @ingroup lang_frame */
const double *seams_frame_positions(const seams_frame *frame);
/** @brief seams_frame_types. @ingroup lang_frame */
const int *seams_frame_types(const seams_frame *frame);
/** @brief seams_frame_ids. @ingroup lang_frame */
const int *seams_frame_ids(const seams_frame *frame);

/** @brief seams_frame_chill_plus. @ingroup lang_order */
int seams_frame_chill_plus(seams_frame *frame, double cutoff, int type,
                           const char *graph, int *labels, int counts9[9]);
/** rule NULL or "" is chill+. "ovito-chill+" is the OVITO 3.17.2 rule.
 * A YAML document or a Lua table fragment with a rule field is accepted.
 * Returns 1 when the rule name is unknown. @ingroup lang_order */
int seams_frame_chill_plus_rule(seams_frame *frame, double cutoff, int type,
                                const char *graph, const char *rule,
                                int *labels, int counts9[9]);
/** @brief seams_frame_chill. @ingroup lang_order */
int seams_frame_chill(seams_frame *frame, double cutoff, int type,
                      const char *graph, int *labels, int counts9[9]);

/** counts[i] is the number of primitive rings of size i.
 * counts_len must be at least max_depth + 1. */
/** @ingroup lang_rings */
int seams_frame_rings(seams_frame *frame, double cutoff, int type,
                      const char *graph, int max_depth, int *counts,
                      int counts_len);

/** Primitive rings of caller CSR. Directed edges are symmetrized. histogram
 * has max_depth + 1 entries. ring_offsets has n_rings + 1 entries and
 * vertices holds their concatenated zero-based indices. NULL row buffers
 * query n_rings and n_vertices; otherwise both capacities are checked.
 * This is the entry used by Python rings_from_csr. @ingroup lang_rings */
int seams_graph_rings(const int *offsets, int n_offsets, const int *columns,
                      int n_columns, int max_depth, int *histogram,
                      int histogram_len, int *ring_offsets,
                      int ring_offsets_cap, int *vertices, int vertices_cap,
                      int *n_rings, int *n_vertices);

/** Signed periodic winding along axis and whole-cycle contractibility for
 * packed rings returned by seams_graph_rings. @ingroup lang_rings */
int seams_frame_ring_winding(const seams_frame *frame,
                             const int *ring_offsets, int n_rings,
                             const int *vertices, int n_vertices, int axis,
                             int *axis_winding, int *contractible);

/** Chordless cycles of one size, from 3 through 12. A cycle is chordless
 * when no bond joins two of its atoms that are not next to each other on
 * it. graph is "cutoff" (NULL or "") or "mutual", the same graphs the C
 * API's neighbour list builds: type > 0 bonds that type, and type <= 0
 * bonds every atom of a frame with more than one species.
 * Each cycle is size atom indices, 0-based over the whole frame. It
 * starts at its lowest index and runs towards the smaller of that atom's
 * two cycle neighbours, and the cycles are sorted, so the rows equal the
 * C API's chordless rows on the same frame. rings may be NULL to count
 * only; otherwise n_rings is required and rings_len is its length in
 * ints. per_atom is NULL or seams_frame_nop ints, the number of cycles
 * through each atom. After seams_frame_select the cycles use the
 * selected atoms only. A short rings buffer returns 1 after n_rings and
 * per_atom are written. @ingroup lang_frame */
int seams_frame_chordless(seams_frame *frame, double cutoff, int type,
                          const char *graph, int size, int *rings,
                          int rings_len, int *n_rings, int *per_atom);

/** Replace coordinates of a frame built by seams_frame_from_arrays.
 * n must be the atom count already stored. The box is nine LAMMPS numbers. @ingroup lang_frame */
int seams_frame_set_positions(seams_frame *frame, const double *xyz, int n,
                             const double *box9);

/** Phase times are seconds. reused is 1 when the wide neighbour list was
 * kept from an earlier call. A null times pointer skips the record. */
typedef struct seams_phase_times {
  double neighbours;
  double chill;
  double rings;
  int reused;
} seams_phase_times;

/** CHILL+ and primitive rings on one neighbour list.
 * skin > 0 builds that list at cutoff+skin and keeps it while every atom's
 * minimum-image displacement stays within skin/2. Bonds are the pairs
 * still inside cutoff, so the labels match a fresh cutoff list.
 * skin <= 0 builds a fresh cutoff list on every call. */
/** @ingroup lang_rings */
int seams_frame_chill_rings(seams_frame *frame, double cutoff, int type,
                           const char *graph, int max_ring, double skin,
                           int *labels, int counts9[9], int *ring_counts,
                           int ring_counts_len, seams_phase_times *times);

/** Build the cutoff neighbour list and discard it. graph NULL means cutoff.
 * seams_frame_neighbour_pairs then returns the directed pairs of that list. @ingroup lang_order */
int seams_frame_neighbours(seams_frame *frame, double cutoff, int type,
                           const char *graph);
/** @brief seams_frame_neighbour_pairs. @ingroup lang_order */
int seams_frame_neighbour_pairs(const seams_frame *frame);

/** H-bond edge export. Rows are sorted by (donor, acceptor, hydrogen).
 * donor and acceptor are frame indices. hydrogen is UINT32_MAX when
 * SEAMS_HBOND_EDGE_WITH_H is not requested. fields is the same bitmask as
 * seams_hbond_edges in seams.h. Any non-null output column needs cap rows.
 * n_out receives the edge count. */
int seams_frame_hbond_edges(seams_frame *frame, double cutoff, int oxygen_type,
                            int hydrogen_type, double dist_cutoff,
                            double angle_cutoff, int fields, uint32_t *donor,
                            uint32_t *acceptor, uint32_t *hydrogen,
                            double *d_DA, double *angle, int cap, int *n_out);

/** Local Steinhardt q_l on the cutoff shell. A null qlbar skips the
 * Lechner-Dellago average; the local values do not use that pass.
 * An empty shell is 0: the harmonic row stays zero and is not divided
 * by the neighbour count. freud and pyscal3 store NaN for that atom. */
/** @ingroup lang_order */
int seams_frame_steinhardt(seams_frame *frame, double cutoff, int type,
                           int order_l, double *ql, double *qlbar);
/** graph is "mutual" (four nearest, the historical default) or "cutoff"
 * (the orthoFullCsr list, the same neighbours as `seams steinhardt`). @ingroup lang_order */
int seams_frame_steinhardt_graph(seams_frame *frame, double cutoff, int type,
                                 int order_l, const char *graph, double *ql,
                                 double *qlbar);

/** Normalized Steinhardt hat{w}_l for l = 4 or 6 on the cutoff shell.
 * Returns 1 on failure. @ingroup lang_order */
int seams_frame_wl(seams_frame *frame, double cutoff, int type, int order_l,
                   double *wl);

/** Local q4 and q6 on the cutoff shell. Degrees 4 and 6 use the
 * solid-harmonic kernel. q4bar and q6bar may be null, and a null
 * pointer skips that Lechner-Dellago average. An empty shell is 0
 * for both the local value and the average. Returns 1 on failure. */
/** @ingroup lang_order */
int seams_frame_q4_q6(seams_frame *frame, double cutoff, int type, double *q4,
                      double *q6, double *q4bar, double *q6bar);

/** One cutoff list. Local q4, local q6, Lechner-Dellago q4-bar and
 * q6-bar, and published CHILL+ labels. q4bar and q6bar may be null.
 * counts9 and labels may be null. Returns 1 on failure. @ingroup lang_order */
int seams_frame_q4_q6_chill(seams_frame *frame, double cutoff, int type,
                            double *q4, double *q6, double *q4bar,
                            double *q6bar, int *labels, int counts9[9]);

/** Errington-Debenedetti q for arbitrary zero-based centre and partner
 * selections. Empty selections mean all atoms. convention 0 is four nearest;
 * convention 1 is every partner within cutoff. Fewer than four writes NaN.
 * @ingroup lang_order */
int seams_frame_tetrahedral(seams_frame *frame, const int *centres,
                            int n_centres, const int *partners,
                            int n_partners, int convention, double cutoff,
                            double *q);

/** Local Structure Index for the same arbitrary selections. Distances inside
 * cutoff are sorted and adjacent radial-gap variance is returned. Fewer than
 * two neighbours writes NaN. @ingroup lang_order */
int seams_frame_lsi(seams_frame *frame, const int *centres, int n_centres,
                    const int *partners, int n_partners, double cutoff,
                    double *lsi);

#ifndef DOXYGEN_SHOULD_SKIP_THIS
/** Repeated prototypes. The documented declarations are above. @ingroup lang_order */
int seams_frame_wl(seams_frame *frame, double cutoff, int type, int order_l,
                   double *wl);
/** @brief seams_frame_q4_q6_chill. @ingroup lang_order */
int seams_frame_q4_q6_chill(seams_frame *frame, double cutoff, int type,
                            double *q4, double *q6, double *q4bar,
                            double *q6bar, int *labels, int counts9[9]);
#endif

/** @brief seams_frame_cages. @ingroup lang_cages */
int seams_frame_cages(seams_frame *frame, double cutoff, int type, int k,
                      const char *graph, int complete, int *per_atom,
                      int *n_hex, int *n_cubic, int *n_water, int *n_both);

/** Ring-size census. spec is "512" or "5:12,6:2". graph is not seeded.
 * The default face-edge window is 2.65,3.15 in the file's coordinate
 * units. seams_frame_signature_edge takes the same string as
 * --cage-edge: NULL or "" keeps that default, and "none" leaves the
 * cage bound to the face census. */
/** @ingroup lang_cages */
int seams_frame_signature(seams_frame *frame, double cutoff, int type,
                          const char *graph, const char *spec, int *n_cages,
                          int *n_atoms);
/** @brief seams_frame_signature_edge. @ingroup lang_cages */
int seams_frame_signature_edge(seams_frame *frame, double cutoff, int type,
                               const char *graph, const char *spec,
                               const char *edge, int *n_cages, int *n_atoms);
/** Atom indices of one matched cage, sorted. vertices may be NULL to
 * query the length. cap is the buffer length. @ingroup lang_cages */
int seams_frame_signature_vertices(seams_frame *frame, double cutoff, int type,
                                   const char *graph, const char *spec,
                                   int cage_index, int *vertices, int cap,
                                   int *n_vertices);
/** Same indices with the face-edge window of seams_frame_signature_edge.
 * NULL or "" keeps 2.65,3.15. "none" leaves the window off. @ingroup lang_cages */
int seams_frame_signature_vertices_edge(seams_frame *frame, double cutoff, int type,
                                        const char *graph, const char *spec,
                                        const char *edge, int cage_index, int *vertices,
                                        int cap, int *n_vertices);

/** @brief seams_frame_rdf_nbins. @ingroup lang_more */
int seams_frame_rdf_nbins(double rmax, int nbins);
/** @brief seams_frame_rdf. @ingroup lang_more */
int seams_frame_rdf(const seams_frame *frame, double rmax, int nbins,
                    int type_i, int type_j, double *r, double *g, int *count,
                    double *volume);

/** @brief seams_frame_rdf2d_nbin. @ingroup lang_more */
int seams_frame_rdf2d_nbin(double cutoff, double binwidth);
/** output_dir NULL skips the topoMonolayer/rdf.dat writer. @ingroup lang_more */
int seams_frame_rdf2d(seams_frame *frame, double cutoff, double binwidth,
                      const char *output_dir, double *r, double *g);

/** @brief seams_prism_row. */
typedef struct seams_prism_row {
  int size;
  int perfect;
  int imperfect;
} seams_prism_row;

/** rows NULL asks for the row count in *n_rows (max_depth - 2). @ingroup lang_more */
int seams_frame_prisms(seams_frame *frame, double cutoff, int type,
                       int max_depth, const char *output_dir,
                       seams_prism_row *rows, int *n_rows);

/** @brief seams_frame_monolayer. @ingroup lang_more */
int seams_frame_monolayer(seams_frame *frame, double cutoff, int type,
                          int max_depth, double sheet_area,
                          const char *output_dir, int *counts, int counts_len);

/** Keep water oxygens for later neighbours, CHILL, CHILL+, rings, cages
 * and Steinhardt. type > 0 is that LAMMPS type. nmol > 0 keeps those
 * molecule ids as well. type == 0 and nmol == 0 infers the oxygen of
 * each multi-atom molecule. A missing type or molecule, or a selection
 * that keeps nothing, returns 1 and sets seams_last_error. Atoms that
 * are not selected stay in the frame and are left unclassified. */
/** @ingroup lang_frame */
int seams_frame_select(seams_frame *frame, int type, const int *mols, int nmol);

/** Same as seams_frame_from_arrays, plus molecule ids. mols may be NULL. @ingroup lang_frame */
seams_frame *seams_frame_from_species(const double *xyz, int n,
                                     const double *box9, const int *types,
                                     const int *ids, const int *mols);

/** Heavy atoms, then a separate hydrogen cloud. Heavies occupy
 * [0, n_heavy). Hydrogens follow in the given order. A null id pointer
 * numbers that cloud from 1 across the frame. A null mol pointer copies
 * the atom id into molID. n_hydrogen may be 0. Returns NULL on failure.
 * @ingroup lang_frame */
seams_frame *seams_frame_from_hydrogen_cloud(
    const double *heavy_xyz, int n_heavy, const int *heavy_types,
    const int *heavy_ids, const int *heavy_mols, const double *hydrogen_xyz,
    int n_hydrogen, const int *hydrogen_types, const int *hydrogen_ids,
    const int *hydrogen_mols, const double *box9);

/** Directed hydrogen bonds. Rows are sorted by donor, hydrogen, acceptor.
 * The angle is acceptor-first: 0 degrees is linear. The acceptor-hydrogen
 * distance is strictly below dist_cutoff. The angle is kept at or below
 * angle_cutoff. The donor-acceptor pair is strictly inside pair_cutoff.
 * exclude_same_mol is 0 or 1. An empty type set writes zero rows. Any
 * non-null output column needs cap rows. Returns 0 on success.
 * @ingroup lang_frame */
int seams_frame_directed_hbonds(seams_frame *frame, double pair_cutoff,
                                double dist_cutoff, double angle_cutoff,
                                const int *donor_types, int n_donor_types,
                                const int *acceptor_types, int n_acceptor_types,
                                int hydrogen_type, int exclude_same_mol,
                                uint32_t *donor, uint32_t *hydrogen,
                                uint32_t *acceptor, double *d_AH, double *d_DA,
                                double *angle, int cap, int *n_out);

/** Per-atom ice label from the shell kernel. labels[i] is 0 ice-Ih,
 * 1 ice-Ic, 2 stacking-disordered, or -1. shell_rmsd is the scale-free
 * residual ceiling (0.12 is the production default). @ingroup lang_more */
int seams_frame_symmetry(seams_frame *frame, double shell_rmsd, int *labels);

/** tau_90 of the H-H vectors between two frames. Writes quiet NaN when no
 * water rotates by 90 degrees. oxygen_type and hydrogen_type select the
 * molecule. Returns 1 when a frame is missing. */
/** @ingroup lang_more */
int seams_frames_jump_tau90(const seams_frame *first, const seams_frame *second,
                            double dt, int oxygen_type, int hydrogen_type,
                            double *tau);

#ifdef __cplusplus
}
#endif

#endif
