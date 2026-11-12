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

#ifdef __cplusplus
}
#endif

#endif
