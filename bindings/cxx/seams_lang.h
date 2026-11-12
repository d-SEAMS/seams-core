#ifndef SEAMS_LANG_H_
#define SEAMS_LANG_H_

#ifdef __cplusplus
extern "C" {
#endif

/* C ABI shared by the Python, R, and Julia fronts.
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
 * graph is NULL or "") or "mutual", the mutual four-nearest list with
 * candidates inside cutoff + 1.5. Returns 0 on success. */

typedef struct seams_frame seams_frame;

const char *seams_lang_version(void);
const char *seams_last_error(void);

int seams_path_nframes(const char *path);

/* Null when path is null or the file cannot be read.
 * seams_last_error() then names the reason. */
seams_frame *seams_frame_read(const char *path, int frame, int type_filter);
seams_frame *seams_frame_from_arrays(const double *xyz, int n,
                                    const double *box9, const int *types,
                                    const int *ids);
void seams_frame_free(seams_frame *frame);

int seams_frame_nop(const seams_frame *frame);
int seams_frame_index(const seams_frame *frame);
void seams_frame_box(const seams_frame *frame, double box9[9]);
const double *seams_frame_positions(const seams_frame *frame);
const int *seams_frame_types(const seams_frame *frame);
const int *seams_frame_ids(const seams_frame *frame);

int seams_frame_chill_plus(seams_frame *frame, double cutoff, int type,
                           const char *graph, int *labels, int counts9[9]);
/* rule NULL or "" is chill+. "ovito-chill+" is the OVITO 3.17.2 rule.
 * A YAML document or a Lua table fragment with a rule field is accepted.
 * Returns 1 when the rule name is unknown. */
int seams_frame_chill_plus_rule(seams_frame *frame, double cutoff, int type,
                                const char *graph, const char *rule,
                                int *labels, int counts9[9]);
int seams_frame_chill(seams_frame *frame, double cutoff, int type,
                      const char *graph, int *labels, int counts9[9]);

/* counts[i] is the number of primitive rings of size i.
 * counts_len must be at least max_depth + 1. */
int seams_frame_rings(seams_frame *frame, double cutoff, int type,
                      const char *graph, int max_depth, int *counts,
                      int counts_len);

/* Replace coordinates of a frame built by seams_frame_from_arrays.
 * n must be the atom count already stored. The box is nine LAMMPS numbers. */
int seams_frame_set_positions(seams_frame *frame, const double *xyz, int n,
                             const double *box9);

/* Phase times are seconds. reused is 1 when the wide neighbour list was
 * kept from an earlier call. A null times pointer skips the record. */
typedef struct seams_phase_times {
  double neighbours;
  double chill;
  double rings;
  int reused;
} seams_phase_times;

/* CHILL+ and primitive rings on one neighbour list.
 * skin > 0 builds that list at cutoff+skin and keeps it while every atom's
 * minimum-image displacement stays within skin/2. Bonds are the pairs
 * still inside cutoff, so the labels match a fresh cutoff list.
 * skin <= 0 builds a fresh cutoff list on every call. */
int seams_frame_chill_rings(seams_frame *frame, double cutoff, int type,
                           const char *graph, int max_ring, double skin,
                           int *labels, int counts9[9], int *ring_counts,
                           int ring_counts_len, seams_phase_times *times);

/* Build the cutoff neighbour list and discard it. graph NULL means cutoff. */
int seams_frame_neighbours(seams_frame *frame, double cutoff, int type,
                           const char *graph);

int seams_frame_steinhardt(seams_frame *frame, double cutoff, int type,
                           int order_l, double *ql, double *qlbar);

int seams_frame_cages(seams_frame *frame, double cutoff, int type, int k,
                      const char *graph, int complete, int *per_atom,
                      int *n_hex, int *n_cubic, int *n_water, int *n_both);

/* Ring-size census. spec is "512" or "5:12,6:2". graph is not seeded.
 * The default face-edge window is 2.65,3.15 in the file's coordinate
 * units. seams_frame_signature_edge takes the same string as
 * --cage-edge: NULL or "" keeps that default, and "none" leaves the
 * cage bound to the face census. */
int seams_frame_signature(seams_frame *frame, double cutoff, int type,
                          const char *graph, const char *spec, int *n_cages,
                          int *n_atoms);
int seams_frame_signature_edge(seams_frame *frame, double cutoff, int type,
                               const char *graph, const char *spec,
                               const char *edge, int *n_cages, int *n_atoms);

int seams_frame_rdf_nbins(double rmax, int nbins);
int seams_frame_rdf(const seams_frame *frame, double rmax, int nbins,
                    int type_i, int type_j, double *r, double *g, int *count,
                    double *volume);

int seams_frame_rdf2d_nbin(double cutoff, double binwidth);
/* output_dir NULL skips the topoMonolayer/rdf.dat writer. */
int seams_frame_rdf2d(seams_frame *frame, double cutoff, double binwidth,
                      const char *output_dir, double *r, double *g);

typedef struct seams_prism_row {
  int size;
  int perfect;
  int imperfect;
} seams_prism_row;

/* rows NULL asks for the row count in *n_rows (max_depth - 2). */
int seams_frame_prisms(seams_frame *frame, double cutoff, int type,
                       int max_depth, const char *output_dir,
                       seams_prism_row *rows, int *n_rows);

int seams_frame_monolayer(seams_frame *frame, double cutoff, int type,
                          int max_depth, double sheet_area,
                          const char *output_dir, int *counts, int counts_len);

/* Keep water oxygens for later neighbours, CHILL, CHILL+, rings, cages
 * and Steinhardt. type > 0 is that LAMMPS type. nmol > 0 keeps those
 * molecule ids as well. type == 0 and nmol == 0 infers the oxygen of
 * each multi-atom molecule. A missing type or molecule, or a selection
 * that keeps nothing, returns 1 and sets seams_last_error. Atoms that
 * are not selected stay in the frame and are left unclassified. */
int seams_frame_select(seams_frame *frame, int type, const int *mols, int nmol);

/* Same as seams_frame_from_arrays, plus molecule ids. mols may be NULL. */
seams_frame *seams_frame_from_species(const double *xyz, int n,
                                     const double *box9, const int *types,
                                     const int *ids, const int *mols);

#ifdef __cplusplus
}
#endif

#endif
