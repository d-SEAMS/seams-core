/* Public C API of Iterative Rotations and Assignments / SOFI.
 * Copyright (c) Miha Gunde and contributors.
 * Apache-2.0. Upstream: mammasmias/IterativeRotationsAssignments.
 * Vendored so a build can link libira without the Fortran module path.
 */
#ifndef IRALIB_H
#define IRALIB_H

#ifdef __cplusplus
extern "C" {
#endif

void libira_cshda(int nat1, int *typ1, double *coords1, int nat2, int *typ2,
                  double *coords2, double thr, int **found, double **dists);

void libira_cshda_from_cost(int n1, int n2, double *cost, int **found,
                            double **dists);

void libira_cshda_pbc(int nat1, int *typ1, double *coords1, int nat2, int *typ2,
                      double *coords2, double *lat2, double thr, int **found,
                      double **dists);

void libira_match(int nat1, int *typ1, double *coords1, int *cand1, int nat2,
                  int *typ2, double *coords2, int *cand2, double km_factor,
                  double **rmat, double **tr, int **perm, double *hd,
                  int *cerr);

void libira_svdrot(int nat1, int *typ1, double *coords1, int nat2, int *typ2,
                   double *coords2, double **rotation, double **translation,
                   int *ierr);

void libira_cshda_svd(int nat1, int *typ1, double *coords1, int nat2, int *typ2,
                      double *coords2, double dthr, int recenter, int **perm,
                      double **dists, double **rmat, double **tr, int *ierr);

void libira_get_version(char *string, long *date);

void libira_compute_all(int nat, int *typ, double *coords, double sym_thr,
                        int prescreen_ih, int *n_mat, double **mat_data,
                        int **perm_data, char **op_data, int **n_data,
                        int **p_data, double **ax_data, double **angle_data,
                        double **dHausdorff_data, char **pg, int *n_prin_ax,
                        double **prin_ax, int *cerr);

void libira_get_symm_ops(int nat, int *typ, double *coords, double sym_thr,
                         int prescreen_ih, int *nmat, double **mat_data,
                         int *cerr);

void libira_get_pg(int n_mat, double *mat_data, char **pg, int *n_prin_ax,
                   double **prin_ax, int verb, int *cerr);

void libira_unique_ax_angle(int nmat, double *mat_data, char **op_data,
                            double **ax_data, double **angle_data, int *cerr);

void libira_analmat(double *mat, char **op, int *n, int *p, double **ax,
                    double *angle, int *cerr);

void libira_get_perm(int nat, int *typ, double *coords, int nmat,
                     double *mat_data, int **perm_data,
                     double **dHausdorff_data);

void libira_get_combos(int nat, int *typ, double *coords, int nmat,
                       double *mat_data, int *nmat_out, double **mat_out,
                       int *cerr);

void libira_try_mat(int nat, int *typ, double *coords, double *rmat, double *dh,
                    int **perm);

void libira_construct_operation(char *op, double *axis, double angle,
                                double **rmat, int *cerr);

void libira_mat_combos(int nmat, double *mat_data, int *nmat_out,
                       double **mat_out);

void libira_ext_bfield(int nmat, double *mat_data, double *b_vec, int *nmat_out,
                       double **mat_out);

void libira_get_err_msg(int ierr, char **msg);

void libira_matrix_distance(const double *mat1, const double *mat2,
                            double *dist);

void libira_check_collinear(int nat, double *coords, int *collinear,
                            double **ax);

#ifdef __cplusplus
}
#endif

#endif
