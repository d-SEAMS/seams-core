#include <stdlib.h>
#include <math.h>
#include <malloc.h>
#include <stdio.h>
#include "externs.h"

void ludcmp(double **a, int n,int *indx,double *d);
void lubksb(double **a,int n, int *indx, double b[]);
void jacobi(matrix a, int n, matrix d, matrix v,int *nrot);

/*
C = A * B  where A is m*n matrix , 
B is k*l matrix, and C is m*l matrix.
C is allocated.  Returns the adress of C .
_______________________________________________________*/

matrix mul_mtx(A,m,n,B,k,l)
int m,n,k,l;
matrix  A, B ;
{
	int i, j , r;
	matrix C, alloc_matrix();
	register double res ;

	if (n!=k ) fprintf(stderr,"Unappropriate dimensions in mul_mtx\n");


 	C = alloc_matrix(m,l) ;


	for (i = 1 ; i <= m ; i++)
		for (j=1 ; j <= l ; j++) {
			res = 0.0 ;
 			for (r=1; r <=n ; r++)
				res += A[i][r]*B[r][j] ;
			C[i][j] = res ;
		}
	return(C) ;
}

/*
C = A + B  where 
A is m*n matrix , 
B is m*n matrix, and C is m*n matrix.
C is allocated.  Returns the adress of C .
_______________________________________________________*/

matrix add_mtx(A,m,n,B,k,l)
int m,n,k,l;
matrix  A, B ;
{
	int i, j;
	matrix C, alloc_matrix();

	if ((m!=k)||(n!=l)) fprintf(stderr,"Unappropriate dimensions in add_mtx\n");


 	C = alloc_matrix(m,n) ;


	for (i = 1 ; i <= m ; i++)
		for (j=1 ; j <= n ; j++) {
			C[i][j] = A[i][j] + B[i][j] ;
		}
	return(C) ;
}


/*
A = A + B  where 
A is m*n matrix , 
B is m*n matrix,
Updates A.
_______________________________________________________*/

void in_place_add_mtx(A,m,n,B,k,l)
int m,n,k,l;
matrix  A, B ;
{
	int i, j;

	if ((m!=k)||(n!=l)) fprintf(stderr,"Unappropriate dimensions in add_mtx\n");


	for (i = 1 ; i <= m ; i++)
		for (j=1 ; j <= n ; j++) {
			A[i][j] +=B[i][j] ;
		}
}

/*
C = A - B  where 
A is m*n matrix , 
B is m*n matrix, and C is m*n matrix.
C is allocated.  Returns the adress of C .
_______________________________________________________*/

matrix sub_mtx(A,m,n,B,k,l)
int m,n,k,l;
matrix  A, B ;
{
        int i, j;
        matrix C, alloc_matrix();

        if ((m!=k)||(n!=l)) fprintf(stderr,"Unappropriate dimensions in add_mtx\n");


        C = alloc_matrix(m,n) ;


        for (i = 1 ; i <= m ; i++)
                for (j=1 ; j <= n ; j++) {
                        C[i][j] = A[i][j] - B[i][j] ;
                }
        return(C) ;
}

/*
C = A * c  where A is m*n matrix , 
c is a double constant.
C is allocated.  Returns the adress of C .
_______________________________________________________*/
matrix mul_const_mtx(A,m,n,c)
int m,n;
double c;
matrix  A;
{
	int i, j ;
	matrix C, alloc_matrix();
 	C = alloc_matrix(m,n) ;
	for (i = 1 ; i <= m ; i++)
	    for (j=1 ; j <= n ; j++)
                C[i][j] = A[i][j]*c;
	return(C) ;
}



/*
A = A*c  where 
A is m*n matrix , 
Updates A.
_______________________________________________________*/

void in_place_mul_const_mtx(A,m,n,c)
int m,n;
double c;
matrix  A;
{
	int i, j;

        for (i = 1 ; i <= m ; i++)
		for (j=1 ; j <= n ; j++) 
			A[i][j] *= c ;
}



/*
C = A x B  where 
A is 3*1 matrix , 
B is 3*1 matrix , 
C is 3*1 matrix , 
C is allocated.  Returns the adress of C .
_______________________________________________________*/
matrix cross_product(A,m,n,B,k,l)
int m,n,k,l;
matrix  A, B ;
{
	matrix C, alloc_matrix();

	if ((n!=1)||(l!=1)) fprintf(stderr,"Unappropriate dimensions in mul_mtx\n");
	if ((m!=3)||(k!=3)) fprintf(stderr,"Unappropriate dimensions in mul_mtx\n");


 	C = alloc_matrix(3,1) ;

        C[1][1] = A[2][1]*B[3][1] - A[3][1]*B[2][1];
        C[2][1] = A[3][1]*B[1][1] - A[1][1]*B[3][1];
        C[3][1] = A[1][1]*B[2][1] - A[2][1]*B[1][1];

	return(C) ;
}

/*
C = <A B>  where 
A is m*1 matrix , 
B is k*1 matrix ,  k=m
C is double.
_______________________________________________________*/
double dot_product(A,m,n,B,k,l)
int m,n,k,l;
matrix  A, B ;
{
	int i;
	double res =0.0;

	if ((n!=1)||(l!=1)) fprintf(stderr,"Unappropriate dimensions in mul_mtx\n");
	if (m!=k) fprintf(stderr,"Unappropriate dimensions in mul_mtx\n");

	for(i=1; i<=m; i++)
            res += A[i][1] * B[i][1];


	return(res) ;
}


/*
C = transpose(A)  where A is m*n matrix . 
 Returns the adress of C .
_______________________________________________________*/

matrix transpose_mtx(A,n,m)
int m,n;
matrix  A;
{
	int i, j ;
	matrix C, alloc_matrix();


 	C = alloc_matrix(m,n) ;


	for (i = 1 ; i <= n ; i++)
		for (j=1 ; j <= m ; j++) 
			C[j][i] = A[i][j] ;
	return(C) ;
}


/***
C = (A)^(-1)
returns the adress of C.
_______________________________________________________***/

matrix inv_mtx(A,n,m)
matrix  A;
int n,m;
{
	int i, j , *indx;
	matrix C, B, alloc_matrix();
	double d ;

	if (m!=n) fprintf(stderr,"in inv_mtx: the matrix is not squared\n") ;

 	C = alloc_matrix(n,n) ;
 	B = alloc_matrix(n,n) ;

	indx = (int *)calloc(n+1,sizeof(int)) ; 

	for (i=1 ; i <=  n ; i++) {
		for (j=1; j <= n ; j++) {
			C[i][j] = 0.0 ;
			B[i][j] = A[i][j] ;
		}
		C[i][i] = 1.0 ;
	}

	ludcmp(B,n,indx,&d);

	for (i=1 ; i <=  n ; i++) 
		lubksb(B,n,indx,C[i]) ;

	
	free(indx);
	free_matrix(B) ;
	B = transpose_mtx(C,n,n) ;
	free_matrix(C);

	return (B) ;
 }


/*
prints the matrix A
_____________________________________________________________*/
void print_mtx(A,m,n)
matrix A ;
int m,n ;
{
	int i, j ;

	for (i=1 ; i <= m  ; i++) {
		for (j=1 ; j <= n  ; j++)
			fprintf(stderr,"  %f",A[i][j]);
		fprintf(stderr,"\n");
	}
}

/*
copy the matrix A
_____________________________________________________________*/
matrix copy_mtx(A,m,n)
matrix A ;
int m,n ;
{
    matrix C, alloc_matrix();

    int i, j ;

    C = alloc_matrix(m,n);
	for (i=1 ; i <= m  ; i++) 
		for (j=1 ; j <= n  ; j++)
			C[i][j] = A[i][j];
	
    return(C);
}

 

/*
alloc m*n array of doubles
_______________________________________________________*/

matrix alloc_matrix(m,n)
int m,n ;
{
	int i ;
	matrix A ;
	double *B ;

	A = (matrix) calloc (m+1,sizeof(B));
	B = (double *) calloc((m+1)*(n+1), sizeof(double)) ;
	for (i=0 ; i<= m ; i++)
			A[i] =  B+i*(n+1);

	return (A);
}

/*
Frees the place of the matrix a.
_______________________________________________________*/

void free_matrix(a)
matrix a;
{
 	free(a[0]) ;
	free(a);
}


/***************************************************************
some utiliyies functions from Numerical-Recepies
****************************************************************/


/**
Given an nXn matrix a, with dimension n, this routie replaces it
by the LU decomposition of a rowwise permutation of itself.
a and n are input. a is output, arranged by superimposed L and U
together. indx is an output vector which records the row permutation 
effected by the partial pivoting. d is output as +-1 depending on
wheather the number of row interchanges was even or odd, resp.
This function is used in combination with LUBKSB to solve linear
equations or invert matrix (see Numerical Recipies pp. 35)
**/

#define TINY 1.0e-20;

void ludcmp(double **a, int n,int *indx,double *d)
{
	int i,imax,j,k;
	double big,dum,sum,temp;
	void nrerror();
	matrix vv, alloc_matrix();

	vv=alloc_matrix(1,n);
	*d=1.0;
	for (i=1;i<=n;i++) {
		big=0.0;
		for (j=1;j<=n;j++)
			if ((temp=fabs(a[i][j])) > big) big=temp;
		if (big == 0.0) nrerror("Singular matrix in routine LUDCMP");
		vv[1][i]=1.0/big;
	}
	for (j=1;j<=n;j++) {
		for (i=1;i<j;i++) {
			sum=a[i][j];
			for (k=1;k<i;k++) sum -= a[i][k]*a[k][j];
			a[i][j]=sum;
		}
		big=0.0;
		for (i=j;i<=n;i++) {
			sum=a[i][j];
			for (k=1;k<j;k++)
				sum -= a[i][k]*a[k][j];
			a[i][j]=sum;
			if ( (dum=vv[1][i]*fabs(sum)) >= big) {
				big=dum;
				imax=i;
			}
		}
		if (j != imax) {
			for (k=1;k<=n;k++) {
				dum=a[imax][k];
				a[imax][k]=a[j][k];
				a[j][k]=dum;
			}
			*d = -(*d);
			vv[1][imax]=vv[1][j];
		}
		indx[j]=imax;
		if (a[j][j] == 0.0) a[j][j]=TINY;
		if (j != n) {
			dum=1.0/(a[j][j]);
			for (i=j+1;i<=n;i++) a[i][j] *= dum;
		}
	}
	free_matrix(vv);
}

#undef TINY




void nrerror(stat)
char *stat;
{
fprintf(stderr,"%s\n",stat);
exit(0);
}

/***
Solves the set of n linear equations AX=B. Here A is input not as the
matrix A but rather as its LU decomposition, determined by the function 
LUDCMP. INDX is input as the permutation vector returned by LUDCMP.
B is input as the right-hand side vector B, and returns with the 
solution vector X. A,N, and INDX are not modified bt this function
and can be left in place for succesive calls with different right-hand
sides B. This function takes into account the possibility that B will
begin with many zero elements, so it is efficient for use in matrix 
inversion
***/

void lubksb(double **a,int n,int *indx, double b[])
{
	int i,ii=0,ip,j;
	double sum;

	for (i=1;i<=n;i++) {
		ip=indx[i];
		sum=b[ip];
		b[ip]=b[i];
		if (ii)
			for (j=ii;j<=i-1;j++) sum -= a[i][j]*b[j];
		else if (sum) ii=i;
		b[i]=sum;
	}
	for (i=n;i>=1;i--) {
		sum=b[i];
		for (j=i+1;j<=n;j++) sum -= a[i][j]*b[j];
		b[i]=sum/a[i][i];
	}
}


/*
This function gets a real symmetric matrix A (nxn) 
and return a matrix (4x4) which is all the eigen vectors (in columns [1-3][i])
and their corresponding eigen values ([4][i]) ordered by size of eigen value.
************************************************************/
matrix all_eigen(A,m,n)
matrix A;
int m,n;
{
        int i,j,nrot,vect1,vect2,vect3;
        double minimum;
        matrix C,V, d, N, alloc_matrix();
;
         
        if (m != n) fprintf(stderr,"ERROR in eigen_val_symm_mtx: not a square matrix.\n");
         
        for (i=1 ; i<=n ; i++)
                for (j=i ; j<=n ; j++)
                        if (A[i][j] != A[j][i])  
                            fprintf(stderr,"ERROR in eigen_val_symm_mtx: non symmetric matrix.\n");

        C = copy_mtx(A,m,n);
        d = alloc_matrix(n,1);
        V = alloc_matrix(n,n);
        N = alloc_matrix(n+1,2*n);

        jacobi(C,n,d,V,&nrot);
        minimum = 10000000.0;
        for(i=1; i<=n; i++)
            if(minimum>d[i][1])
            {   minimum = d[i][1];
                vect1 = i;
            }
        for(i=1; i<=n; i++)
            N[i][1] = V[i][vect1];
        N[n+1][1] = minimum;

        if(d[(vect1%3)+1][1]<d[(((vect1%3)+1)%3)+1][1])
        {   vect2 =   (vect1%3)+1;
            vect3 = (((vect1%3)+1)%3)+1;
        }
        else
        {   vect2 = (((vect1%3)+1)%3)+1;
            vect3 =   (vect1%3)+1;
        }

        for(i=1; i<=n; i++)
        {   N[i][2] = V[i][vect2];
            N[i][3] = V[i][vect3];
        }
        N[n+1][2] = d[vect2][1];
        N[n+1][3] = d[vect3][1];

        for(i=1; i<=n; i++)
        {   N[1][i+n] = -1.0*N[1][i];
            N[2][i+n] = -1.0*N[2][i];
            N[3][i+n] = -1.0*N[3][i];
            N[4][i+n] = -1.0*N[4][i];
        }

        free_matrix(C);
        free_matrix(d);
        free_matrix(V);

        return(N);
}



/*
This function gets a real symmetric real symmetric matrix A (nxn) 
and return a double which is the smallest eigen value and returns in the 
vector n (nx1) the corresponding eigen vector.
************************************************************/
matrix smallest_eigen(A,m,n,g)
matrix A;
int m,n;
double *g;
{
        int i,j,nrot,vect;
        double minimum;
        matrix C,V, d, N, alloc_matrix();
;
         
        if (m != n) fprintf(stderr,"ERROR in eigen_val_symm_mtx: not a square matrix.\n");
         
        for (i=1 ; i<=n ; i++)
                for (j=i ; j<=n ; j++)
                        if (A[i][j] != A[j][i])  
                            fprintf(stderr,"ERROR in eigen_val_symm_mtx: non symmetric matrix.\n");

        N = alloc_matrix(n,1);
        C = copy_mtx(A,m,n);
        d = alloc_matrix(n,1);
        V = alloc_matrix(n,n);

        jacobi(C,n,d,V,&nrot);
        minimum = 10000000.0;
        for(i=1; i<=n; i++)
            if(minimum>d[i][1])
            {   minimum = d[i][1];
                vect = i;
            }
        *g = minimum;
        for(i=1; i<=n; i++)
            N[i][1] = V[i][vect];

        
        free_matrix(C);
        free_matrix(d);
        free_matrix(V);

        return(N);
}



/* 
Computes all eigenvalues and eigenvectors of a real symmetric matrix
a[1..n,1..n]. On output, elements of a above the diagonal are destroyed.
d[1..n] returns the eigenvalues of a. v[1..n,1..n] is a matrix whose
columns contain, on output, the normalized eigenvectors of a.
nrot returns the number of jacobi rotations that where required.
----------------------------------------------------------------------**/

#define ROTATE(a,i,j,k,l) g=a[i][j];h=a[k][l];a[i][j]=g-s*(h+g*tau);\
        a[k][l]=h+s*(g-h*tau);

void jacobi(matrix a, int n, matrix d, matrix v,int *nrot)
{
        int j,iq,ip,i;
        double tresh,theta,tau,t,sm,s,h,g,c;
        matrix b=NULL,z=NULL, alloc_matrix();
        void nrerror();

        b = alloc_matrix(n,1); 
        z = alloc_matrix(n,1); 
        for (ip=1;ip<=n;ip++) {
                for (iq=1;iq<=n;iq++) v[ip][iq]=0.0;
                v[ip][ip]=1.0;
        }
        for (ip=1;ip<=n;ip++) {
                b[ip][1]=d[ip][1]=a[ip][ip];
                z[ip][1]=0.0;
        }
        *nrot=0;
        for (i=1;i<=50;i++) {
                sm=0.0;
                for (ip=1;ip<=n-1;ip++) {
                        for (iq=ip+1;iq<=n;iq++)
                                sm += fabs(a[ip][iq]);
                }
                if (sm == 0.0) {
                        free_matrix(z);
                        free_matrix(b);
                        return;
                }
                if (i < 4)
                        tresh=0.2*sm/(n*n);
                else
                        tresh=0.0;
                for (ip=1;ip<=n-1;ip++) {
                        for (iq=ip+1;iq<=n;iq++) {
                                g=100.0*fabs(a[ip][iq]);
                                if (i > 4 && fabs(d[ip][1])+g == fabs(d[ip][1])
                                        && fabs(d[iq][1])+g == fabs(d[iq][1]))
                                        a[ip][iq]=0.0;
                                else if (fabs(a[ip][iq]) > tresh) {
                                        h=d[iq][1]-d[ip][1];
                                        if (fabs(h)+g == fabs(h))
                                                t=(a[ip][iq])/h;
                                        else {
                                                theta=0.5*h/(a[ip][iq]);
                                                t=1.0/(fabs(theta)+sqrt(1.0+theta*theta));
                                                if (theta < 0.0) t = -t;
                                        }
                                        c=1.0/sqrt(1+t*t);
                                        s=t*c;
                                        tau=s/(1.0+c);
                                        h=t*a[ip][iq];
                                        z[ip][1] -= h;
                                        z[iq][1] += h;
                                        d[ip][1] -= h;
                                        d[iq][1] += h;
                                        a[ip][iq]=0.0;                                        for (j=1;j<=ip-1;j++) {
                                                ROTATE(a,j,ip,j,iq)
                                        }
                                        for (j=ip+1;j<=iq-1;j++) {
                                                ROTATE(a,ip,j,j,iq)
                                        }
                                        for (j=iq+1;j<=n;j++) {
                                                ROTATE(a,ip,j,iq,j)
                                        }
                                        for (j=1;j<=n;j++) {
                                                ROTATE(v,j,ip,j,iq)
                                        }
                                        ++(*nrot);
                                }
                        }
                }
                for (ip=1;ip<=n;ip++) {
                        b[ip][1] += z[ip][1];
                        d[ip][1]=b[ip][1];
                        z[ip][1]=0.0;
                }
        }
        nrerror("Too many iterations in routine JACOBI");
}

#undef ROTATE
