#ifndef SEAMS_C_API_H_
#define SEAMS_C_API_H_

#ifdef __cplusplus
extern "C" {
#endif

/* C ABI for native callers and for an emscripten/WASM build.
 * xyz is n*3, box is ox,oy,oz,lx,ly,lz (ortho). labels[i] is
 * 0 cubic, 1 hexagonal, 2 water, 3 interfacial, 4 clathrate,
 * 5 interClathrate, 6 unclassified, 7 reCubic, 8 reHex.
 * Returns 0 on success.
 * seams_chill_plus bonds every pair closer than 3.5 A (the cutoff list).
 * seams_chill_plus_graph takes the cutoff and graph: "cutoff" (also NULL
 * or "") or "mutual", the mutual four nearest inside cutoff + 1.5. */
int seams_chill_plus(const double *xyz, int n, const double *box, int *labels);
int seams_chill_plus_graph(const double *xyz, int n, const double *box,
                           double cutoff, const char *graph, int *labels);
/* rule is "chill+" (NULL or "") or "ovito-chill+". A YAML or Lua
 * document with a rule field is accepted. Returns 1 on a bad rule. */
int seams_chill_plus_rule(const double *xyz, int n, const double *box,
                          double cutoff, const char *graph, const char *rule,
                          int *labels);

/* box9 is xlo, ylo, zlo, lx, ly, lz, xy, xz, yz. lx, ly, lz are the
 * edge lengths xhi-xlo, not the bound spans, and xlo, ylo, zlo are the
 * edge origin. A zero tilt matches seams_chill_plus_graph. */
int seams_chill_plus_graph_lammps(const double *xyz, int n,
                                  const double *box9, double cutoff,
                                  const char *graph, int *labels);

/* Local symmetry. points is n*3. Template ids:
 * 0 ice-Ih, 1 ice-Ic, 2 ice-Isd, 3 5^12, 4 5^12 6^2, 5 5^12 6^4,
 * 6 ring-5, 7 ring-6.
 * permutation[i] is the environment index matched to template atom i.
 * rotation9 is row-major and, with translation3, takes that environment
 * atom onto the template. Returns 0 on success. */
int seams_symmetry_match(const double *points, int n, int *template_id,
                         double *rmsd, double *hausdorff, double *rotation9,
                         double *translation3, int *permutation);
/* Returns 2 when SOFI is not linked. */
int seams_symmetry_point_group(const double *points, int n, double sym_thr,
                               char *symbol, int symbol_len, int *n_operations);
/* lattice9 is three row vectors (ax,ay,az, bx,by,bz, cx,cy,cz).
 * Returns 2 when spglib is not linked. */
int seams_symmetry_space_group(const double *lattice9, const double *cartesian,
                               int n, double symprec, int *number, char *symbol,
                               int symbol_len);
/* Writes the ideal template. points may be null to query *n_out. */
int seams_symmetry_template(int id, double *points, int cap, int *n_out);
/* which: 0 ice Ih, 1 ice Ic, 2 stacking-fault cell. */
int seams_symmetry_crystal(int which, double *lattice9, double *xyz, int cap,
                           int *n_out);
int seams_symmetry_backend(char *name, int name_len);
/* Per-atom ice label. lattice9 is three row vectors. labels[i] is
 * 0 ice-Ih, 1 ice-Ic, 2 ice-Isd, or -1. The second-shell RMSD ceiling
 * is the scale-free ceiling kShellRmsd (0.12). */
int seams_symmetry_classify(const double *xyz, int n, const double *lattice9,
                            int *labels);
/* Same labels with an explicit scale-free residual ceiling. */
int seams_symmetry_classify_tol(const double *xyz, int n, const double *lattice9,
                                double shell_rmsd, int *labels);

/* Mixed frames. types and mols may be NULL (every atom is then one
 * species, mol id = index + 1). oxygen_type > 0 keeps that LAMMPS type.
 * mol_mode 0 ignores molecules, 1 keeps the water oxygen of each
 * multi-atom molecule, 2 keeps mol_ids[0..n_mol_ids). A mixed frame
 * with mol_mode 0 and oxygen_type 0 is an error. Labels of atoms that
 * are not water oxygens are unclassified (6). ion_ice, ion_liquid and
 * ion_front may be NULL. ion_front is the nearest ice oxygen, or -1.
 * seams_c_error() is the message after a non-zero return.
 * On a single-species cloud the labels match seams_chill_plus_graph_lammps. */
int seams_chill_plus_species(const double *xyz, int n, const double *box9,
                             const int *types, const int *mols,
                             int oxygen_type, int mol_mode, const int *mol_ids,
                             int n_mol_ids, double cutoff, const char *graph,
                             int *labels, int *ion_ice, int *ion_liquid,
                             double *ion_front);

/* degree[i] is the oxygen-oxygen neighbour count inside cutoff, or -1
 * when atom i is not a selected water oxygen. */
int seams_oxygen_cutoff_degree(const double *xyz, int n, const double *box9,
                               const int *types, const int *mols,
                               int oxygen_type, int mol_mode,
                               const int *mol_ids, int n_mol_ids, double cutoff,
                               int *degree);

const char *seams_c_error(void);

/* Local q4 and q6 on the cutoff shell. Degrees 4 and 6 use the
 * solid-harmonic kernel. xyz is n*3. box is ox,oy,oz,lx,ly,lz.
 * q4 and q6 are n doubles. q4bar and q6bar may be null; a null pointer
 * skips that Lechner-Dellago average. Returns 0 on success. */
int seams_q4_q6(const double *xyz, int n, const double *box, double cutoff,
                double *q4, double *q6, double *q4bar, double *q6bar);
/* box9 is xlo,ylo,zlo,lx,ly,lz,xy,xz,yz, the same cell as
 * seams_chill_plus_graph_lammps. */
int seams_q4_q6_lammps(const double *xyz, int n, const double *box9,
                       double cutoff, double *q4, double *q6, double *q4bar,
                       double *q6bar);

/* box9 is xlo, ylo, zlo, lx, ly, lz, xy, xz, yz, the same edges as
 * seams_chill_plus_graph_lammps. The cutoff list is orthoFullCsr: cloud
 * indices, no self entry, both directions, squared distance strictly
 * below cutoff squared.
 * offsets has n + 1 entries. When cols is null, *n_cols receives the
 * pair count and the call returns 0. When cols_cap is short, *n_cols is
 * that count and the call returns 2. */
int seams_cutoff_csr(const double *xyz, int n, const double *box9, double cutoff,
                     int *offsets, int *cols, int cols_cap, int *n_cols);
/* Local Steinhardt q_l on that cutoff list. order_l is 3, 4, 6, 8, or 12. */
int seams_ql_cutoff(const double *xyz, int n, const double *box9, double cutoff,
                    int order_l, double *ql);
/* Diamond common-neighbour labels. 1 cubic, 0 hexagonal, -1 otherwise. */
int seams_diamond_cna(const double *xyz, int n, const double *box9, double cutoff,
                      int *labels);
/* counts[k] is the number of primitive rings of size k.
 * counts_len must be at least max_depth + 1. max_depth is at least 3. */
int seams_primitive_rings(const double *xyz, int n, const double *box9,
                          double cutoff, int max_depth, int *counts, int counts_len);

#ifdef __cplusplus
}
#endif

#endif
