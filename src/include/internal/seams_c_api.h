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

#ifdef __cplusplus
}
#endif

#endif
