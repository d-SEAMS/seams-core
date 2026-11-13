#ifndef EXTERNS_H
#define EXTERNS_H

/* DEFINITIONS */
#define ADJ(xx,yy)    adjacency[xx][yy]
#define ADJ_OLD(xx,yy)    adjacency_old[xx][yy]
#define PRINT_ERR(s)   {fprintf(stderr,s); exit (0);}
#define MAX(nn1,nn2) ((nn1)>=(nn2) ? nn1 : nn2)
#define MAXPOINTS   320

typedef  double  **matrix ;

/* VARIABLES */
extern int count;
extern int adjacency[MAXPOINTS][MAXPOINTS];
extern int adjacency_old[MAXPOINTS][MAXPOINTS];
extern double pdir[2][3];
extern double pdraw[MAXPOINTS][3];
extern double pndraw[MAXPOINTS][3];
extern matrix dir;
extern matrix current_axis;
extern matrix parray[MAXPOINTS];
extern int best_q[2][MAXPOINTS+2];

int svdcmp(matrix a,int m,int n,matrix w,matrix v);

matrix alloc_matrix(int m,int n);
void free_matrix(matrix a);
matrix mul_const_mtx(matrix A,int m, int n,double c);
matrix mul_mtx(matrix A,int m,int n,matrix B,int k,int l);
matrix transpose_mtx(matrix A, int m, int n);
void in_place_add_mtx(matrix A,int m,int n,matrix B,int k,int l);

/* PROCEDURE */

int read_input(int number_case, char *fname, char *shape_file);

double dot_product(matrix A, int m,int n,matrix B, int k,int l);

#endif
