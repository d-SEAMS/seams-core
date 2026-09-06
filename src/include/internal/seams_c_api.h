#ifndef SEAMS_C_API_H_
#define SEAMS_C_API_H_

#ifdef __cplusplus
extern "C" {
#endif

/* C ABI for native callers and for an emscripten/WASM build.
 * xyz is n*3, box is ox,oy,oz,lx,ly,lz (ortho). labels[i] is
 * 0 cubic, 1 hexagonal, 2 water, 3 interfacial, 4 clathrate,
 * 5 interClathrate, 6 unclassified, 7 reCubic, 8 reHex.
 * Returns 0 on success. */
int seams_chill_plus(const double *xyz, int n, const double *box, int *labels);

#ifdef __cplusplus
}
#endif

#endif
