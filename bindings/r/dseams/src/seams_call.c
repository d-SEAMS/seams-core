/* Copyright (c) d-SEAMS developers */
/* .Call shim over include/seams.h. No C++ and no engine headers. */

#include <stdio.h>

#include <R.h>
#include <Rinternals.h>

#include "seams.h"

static seams_system *sys_of(SEXP ext) {
  if (TYPEOF(ext) != EXTPTRSXP) {
    error("system handle");
  }
  seams_system *sys = (seams_system *)R_ExternalPtrAddr(ext);
  if (sys == NULL) {
    error("empty system");
  }
  return sys;
}

static void sys_fin(SEXP ext) {
  seams_system *sys = (seams_system *)R_ExternalPtrAddr(ext);
  if (sys != NULL) {
    seams_system_destroy(sys);
  }
  R_ClearExternalPtr(ext);
}

static void die(const void *handle) {
  const char *text = seams_handle_error(handle);
  error("%s", text == NULL ? "" : text);
}

static void copy_error(const void *handle, char *saved, size_t cap) {
  const char *text = seams_handle_error(handle);
  snprintf(saved, cap, "%s", text == NULL ? "" : text);
}

static SEXP real_n(R_xlen_t n) { return allocVector(REALSXP, n); }

SEXP seams_r_system(void) {
  seams_system *sys = seams_system_create();
  if (sys == NULL) {
    error("null system");
  }
  SEXP ext = PROTECT(R_MakeExternalPtr(sys, R_NilValue, R_NilValue));
  R_RegisterCFinalizerEx(ext, sys_fin, TRUE);
  UNPROTECT(1);
  return ext;
}

SEXP seams_r_read(SEXP ext, SEXP path, SEXP frame, SEXP type) {
  seams_system *sys = sys_of(ext);
  const char *text = Rf_translateCharUTF8(STRING_ELT(path, 0));
  if (seams_system_read(sys, text, Rf_asInteger(frame), Rf_asInteger(type)) != SEAMS_OK) {
    die(sys);
  }
  return ScalarInteger(seams_system_count(sys));
}

SEXP seams_r_count(SEXP ext) { return ScalarInteger(seams_system_count(sys_of(ext))); }

SEXP seams_r_threads(SEXP n) {
  int threads = Rf_asInteger(n);
  if (seams_set_num_threads(threads) != SEAMS_OK) {
    error("%s", seams_last_error());
  }
  return ScalarInteger(seams_get_num_threads());
}

SEXP seams_r_chill(SEXP ext, SEXP cutoff, SEXP type, SEXP graph, SEXP rule) {
  seams_system *sys = sys_of(ext);
  const int n = seams_system_count(sys);
  SEXP labels = PROTECT(allocVector(INTSXP, n));
  SEXP counts = PROTECT(allocVector(INTSXP, 9));
  if (seams_chill(sys, Rf_asReal(cutoff), Rf_asInteger(type), Rf_asInteger(graph), Rf_asInteger(rule),
                  INTEGER(labels), (size_t)n, INTEGER(counts)) != SEAMS_OK) {
    die(sys);
  }
  UNPROTECT(2);
  return counts;
}

SEXP seams_r_steinhardt(SEXP ext, SEXP cutoff, SEXP order, SEXP average) {
  seams_system *sys = sys_of(ext);
  const int n = seams_system_count(sys);
  SEXP ql = PROTECT(real_n(n));
  SEXP bar = PROTECT(real_n(Rf_asInteger(average) ? n : 0));
  if (seams_steinhardt(sys, Rf_asReal(cutoff), 0, Rf_asInteger(order), REAL(ql), (size_t)n,
                       Rf_asInteger(average) ? REAL(bar) : NULL, (size_t)n) != SEAMS_OK) {
    die(sys);
  }
  SEXP out = PROTECT(allocVector(VECSXP, 2));
  SET_VECTOR_ELT(out, 0, ql);
  SET_VECTOR_ELT(out, 1, bar);
  UNPROTECT(3);
  return out;
}

SEXP seams_r_cages(SEXP ext, SEXP cutoff, SEXP graph, SEXP k, SEXP complete) {
  seams_system *sys = sys_of(ext);
  seams_cage_set *cages = NULL;
  if (seams_cages(sys, Rf_asReal(cutoff), 0, Rf_asInteger(graph), Rf_asInteger(k),
                  Rf_asInteger(complete), &cages) != SEAMS_OK) {
    die(sys);
  }
  int n_hex = 0, n_cubic = 0, n_water = 0, n_both = 0;
  const int code = seams_cages_counts(cages, &n_hex, &n_cubic, &n_water, &n_both);
  char saved[512];
  if (code != SEAMS_OK) {
    copy_error(cages, saved, sizeof(saved));
  }
  seams_cages_destroy(cages);
  if (code != SEAMS_OK) {
    error("%s", saved);
  }
  SEXP out = PROTECT(allocVector(INTSXP, 4));
  INTEGER(out)[0] = n_hex;
  INTEGER(out)[1] = n_cubic;
  INTEGER(out)[2] = n_water;
  INTEGER(out)[3] = n_both;
  UNPROTECT(1);
  return out;
}

SEXP seams_r_signature(SEXP ext, SEXP cutoff, SEXP graph, SEXP spec) {
  seams_system *sys = sys_of(ext);
  size_t n_cages = 0, n_atoms = 0;
  if (seams_signature_cages(sys, Rf_asReal(cutoff), 0, Rf_asInteger(graph),
                            Rf_translateCharUTF8(STRING_ELT(spec, 0)), &n_cages, &n_atoms) != SEAMS_OK) {
    die(sys);
  }
  SEXP out = PROTECT(allocVector(INTSXP, 2));
  INTEGER(out)[0] = (int)n_cages;
  INTEGER(out)[1] = (int)n_atoms;
  UNPROTECT(1);
  return out;
}

SEXP seams_r_knn(SEXP ext, SEXP k, SEXP cutoff, SEXP mutual) {
  seams_system *sys = sys_of(ext);
  seams_neighbours *nb = NULL;
  if (seams_neighbours_knn(sys, Rf_asInteger(k), Rf_asReal(cutoff), 0, Rf_asInteger(mutual), &nb) !=
      SEAMS_OK) {
    die(sys);
  }
  size_t n_offsets = 0, nnz = 0;
  const int code = seams_neighbours_copy(nb, NULL, 0, &n_offsets, NULL, 0, &nnz);
  char saved[512];
  if (code != SEAMS_OK) {
    copy_error(nb, saved, sizeof(saved));
  }
  seams_neighbours_destroy(nb);
  if (code != SEAMS_OK) {
    error("%s", saved);
  }
  return ScalarInteger((int)nnz);
}

SEXP seams_r_rdf(SEXP ext, SEXP rmax, SEXP bins) {
  seams_system *sys = sys_of(ext);
  size_t n_out = 0;
  double volume = 0.0;
  if (seams_rdf(sys, Rf_asReal(rmax), Rf_asInteger(bins), 0, 0, NULL, NULL, NULL, 0, &n_out,
                &volume) != SEAMS_OK) {
    die(sys);
  }
  SEXP radius = PROTECT(real_n((R_xlen_t)n_out));
  SEXP g = PROTECT(real_n((R_xlen_t)n_out));
  if (seams_rdf(sys, Rf_asReal(rmax), Rf_asInteger(bins), 0, 0, REAL(radius), REAL(g), NULL, n_out,
                &n_out, NULL) != SEAMS_OK) {
    die(sys);
  }
  SEXP vol = PROTECT(ScalarReal(volume));
  SEXP out = PROTECT(allocVector(VECSXP, 3));
  SET_VECTOR_ELT(out, 0, radius);
  SET_VECTOR_ELT(out, 1, g);
  SET_VECTOR_ELT(out, 2, vol);
  UNPROTECT(4);
  return out;
}

SEXP seams_r_running_cn(SEXP ext, SEXP rmax, SEXP bins) {
  seams_system *sys = sys_of(ext);
  size_t n_out = 0;
  if (seams_running_cn(sys, Rf_asReal(rmax), Rf_asInteger(bins), 0, 0, NULL, NULL, 0, &n_out) !=
      SEAMS_OK) {
    die(sys);
  }
  SEXP radius = PROTECT(real_n((R_xlen_t)n_out));
  SEXP cn = PROTECT(real_n((R_xlen_t)n_out));
  if (seams_running_cn(sys, Rf_asReal(rmax), Rf_asInteger(bins), 0, 0, REAL(radius), REAL(cn), n_out,
                       &n_out) != SEAMS_OK) {
    die(sys);
  }
  UNPROTECT(2);
  return cn;
}

SEXP seams_r_coordination(SEXP ext, SEXP rmax, SEXP bins) {
  seams_system *sys = sys_of(ext);
  double cn = 0.0;
  if (seams_coordination(sys, Rf_asReal(rmax), Rf_asInteger(bins), 0, 0, &cn) != SEAMS_OK) {
    die(sys);
  }
  return ScalarReal(cn);
}

SEXP seams_r_hbonds(SEXP ext, SEXP cutoff, SEXP oxygen, SEXP hydrogen, SEXP dist, SEXP angle) {
  seams_system *sys = sys_of(ext);
  int n_bonds = 0;
  if (seams_hbonds(sys, Rf_asReal(cutoff), Rf_asInteger(oxygen), Rf_asInteger(hydrogen),
                   Rf_asReal(dist), Rf_asReal(angle), &n_bonds) != SEAMS_OK) {
    die(sys);
  }
  return ScalarInteger(n_bonds);
}

SEXP seams_r_density(SEXP ext, SEXP type, SEXP bins, SEXP axis) {
  seams_system *sys = sys_of(ext);
  size_t n_out = 0;
  if (seams_density_axis(sys, Rf_asInteger(type), Rf_asInteger(bins), Rf_asInteger(axis), NULL, NULL,
                         0, &n_out) != SEAMS_OK) {
    die(sys);
  }
  SEXP z = PROTECT(real_n((R_xlen_t)n_out));
  SEXP rho = PROTECT(real_n((R_xlen_t)n_out));
  if (seams_density_axis(sys, Rf_asInteger(type), Rf_asInteger(bins), Rf_asInteger(axis), REAL(z),
                         REAL(rho), n_out, &n_out) != SEAMS_OK) {
    die(sys);
  }
  UNPROTECT(2);
  return rho;
}

SEXP seams_r_fingerprint(SEXP ext, SEXP cutoff, SEXP hops) {
  seams_system *sys = sys_of(ext);
  seams_topo *fp = NULL;
  if (seams_fingerprint(sys, Rf_asReal(cutoff), 0, 4, SEAMS_GRAPH_CUTOFF, Rf_asInteger(hops), 0,
                        &fp) != SEAMS_OK) {
    die(sys);
  }
  size_t n_out = 0;
  int counts[13] = {0};
  int code = seams_fingerprint_key(fp, NULL, 0, &n_out);
  if (code == SEAMS_OK) {
    code = seams_fingerprint_rings(fp, counts, 13);
  }
  char saved[512];
  if (code != SEAMS_OK) {
    copy_error(fp, saved, sizeof(saved));
  }
  seams_fingerprint_destroy(fp);
  if (code != SEAMS_OK) {
    error("%s", saved);
  }
  int ring_sum = 0;
  for (int i = 0; i < 13; ++i) {
    ring_sum += counts[i];
  }
  SEXP out = PROTECT(allocVector(INTSXP, 2));
  INTEGER(out)[0] = (int)n_out;
  INTEGER(out)[1] = ring_sum;
  UNPROTECT(1);
  return out;
}

SEXP seams_r_pairs(SEXP ext, SEXP spec) {
  seams_system *sys = sys_of(ext);
  size_t n_pairs = 0;
  int n_cation = 0, n_anion = 0;
  if (seams_contact_pairs(sys, Rf_translateCharUTF8(STRING_ELT(spec, 0)), &n_pairs, &n_cation,
                          &n_anion) != SEAMS_OK) {
    die(sys);
  }
  SEXP out = PROTECT(allocVector(INTSXP, 3));
  INTEGER(out)[0] = (int)n_pairs;
  INTEGER(out)[1] = n_cation;
  INTEGER(out)[2] = n_anion;
  UNPROTECT(1);
  return out;
}

SEXP seams_r_domains(SEXP ext, SEXP cutoff, SEXP spec, SEXP subset) {
  seams_system *sys = sys_of(ext);
  int n_subset = 0, largest = 0;
  double percolation = 0.0;
  if (seams_domains(sys, Rf_asReal(cutoff), Rf_translateCharUTF8(STRING_ELT(spec, 0)),
                    Rf_asInteger(subset), &n_subset, &largest, &percolation) != SEAMS_OK) {
    die(sys);
  }
  SEXP out = PROTECT(real_n(3));
  REAL(out)[0] = n_subset;
  REAL(out)[1] = largest;
  REAL(out)[2] = percolation;
  UNPROTECT(1);
  return out;
}

SEXP seams_r_prisms(SEXP ext, SEXP cutoff) {
  seams_system *sys = sys_of(ext);
  SEXP perfect = PROTECT(allocVector(INTSXP, 7));
  SEXP imperfect = PROTECT(allocVector(INTSXP, 7));
  if (seams_prisms(sys, Rf_asReal(cutoff), 0, INTEGER(perfect), INTEGER(imperfect), 7) != SEAMS_OK) {
    die(sys);
  }
  UNPROTECT(2);
  return perfect;
}

SEXP seams_r_monolayer(SEXP ext, SEXP cutoff, SEXP depth, SEXP area) {
  seams_system *sys = sys_of(ext);
  const int max_depth = Rf_asInteger(depth);
  SEXP counts = PROTECT(allocVector(INTSXP, max_depth + 1));
  if (seams_monolayer(sys, Rf_asReal(cutoff), 0, max_depth, Rf_asReal(area), INTEGER(counts),
                      (size_t)(max_depth + 1)) != SEAMS_OK) {
    die(sys);
  }
  UNPROTECT(1);
  return counts;
}

SEXP seams_r_vertices(SEXP ext, SEXP cutoff, SEXP depth, SEXP index) {
  seams_system *sys = sys_of(ext);
  seams_rings *rings = NULL;
  if (seams_rings_find(sys, Rf_asReal(cutoff), 0, SEAMS_GRAPH_CUTOFF, Rf_asInteger(depth), &rings) !=
      SEAMS_OK) {
    die(sys);
  }
  size_t n_out = 0;
  int code = seams_rings_vertices(rings, Rf_asInteger(index), NULL, 0, &n_out);
  SEXP atoms = PROTECT(allocVector(INTSXP, (R_xlen_t)n_out));
  if (code == SEAMS_OK) {
    code = seams_rings_vertices(rings, Rf_asInteger(index), INTEGER(atoms), n_out, &n_out);
  }
  char saved[512];
  if (code != SEAMS_OK) {
    copy_error(rings, saved, sizeof(saved));
  }
  seams_rings_destroy(rings);
  if (code != SEAMS_OK) {
    error("%s", saved);
  }
  UNPROTECT(1);
  return atoms;
}

SEXP seams_r_site(SEXP ext, SEXP spec) {
  seams_system *sys = sys_of(ext);
  const int n = seams_system_count(sys);
  SEXP kinds = PROTECT(allocVector(INTSXP, n));
  if (seams_site_table(sys, Rf_translateCharUTF8(STRING_ELT(spec, 0)), INTEGER(kinds), (size_t)n) !=
      SEAMS_OK) {
    die(sys);
  }
  UNPROTECT(1);
  return kinds;
}
