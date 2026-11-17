extern "C" {
#include <math.h>
#include "externs.h"
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include "externs.h"
#include "Molecule.h"
#include "permuter.h"
}

#include <openbabel/mol.h>
#include "babelAdapter.h"

#define MAX_STRING_LEN 1024
#define MAX_NUM_ARG 64
#define ARG_LENGTH 64
#define BUF_LEN 64

#define MAXDOUBLE   100000000.0
#define TWOPI       6.283185307
#define HALFPI      1.570796327
#define PI          3.141592654
#define SMALL       0.0000001
#define SCALE       1.0
#define FALSE       0
#define TRUE        1
#define ABS(x)      ((x) > 0.0 ? (x) : -1.0*(x))
#define SQR(x)      ((x) * (x))
#define irint(xx)   (int)rint(xx)

#define SHAPE	    0
#define  B4         1
#define AB4         2
#define  B5         3
#define AB5         4
#define  B6         5
#define AB6         6
#define  B8         7
#define AB8         8
#define  B12        9
#define AB12       10
#define  B20       11
#define AB20       12
#define  B60       13
#define AB60       14
#define  SQ5       15
#define  B5_       16
#define AB5_       17
#define  PR        19
#define APR        20
#define  PR_EQ     21
#define APR_EQ     22
#define  D3H       23
#define AD3H       24


double x_avg      = 0.0,
       y_avg      = 0.0,
       z_avg      = 0.0;
double s          = MAXDOUBLE;
double shift      = 0.0;
double norm       = 1.0;
double shift_symm = 0.0;
double rms_symm   = 0.0;
double edge;
double norma      = 1.0;
int count         = 0;

int valency[MAXPOINTS];
int q[2][MAXPOINTS+2];
int match[MAXPOINTS];
double current_s;
int best_permute[MAXPOINTS];
int best_perm[MAXPOINTS];
int best_q[2][MAXPOINTS+2];
matrix parray[MAXPOINTS];
matrix cbarray[MAXPOINTS];
OBMol inMol;
Molecule* mol;
OBMol inShape;
Molecule* shape; 

int **permutations = NULL;
int numPermutations = 0;
char *permsFile = NULL;
int allowImproperRotations = 0;



int adj[MAXPOINTS][MAXPOINTS];
int adjacency[MAXPOINTS][MAXPOINTS];
int adjacency_old[MAXPOINTS][MAXPOINTS];

int tail = -1;
matrix dir;
matrix rotation;

double pdraw[MAXPOINTS][3];
double pndraw[MAXPOINTS][3];
char file_buffer[2000][256];
int file_n_lines=0;
char file_out[1024];
char file_shape_out[1024];
char file_shape_in[1024];
char file_model[1024];
void eval_cube_symm();
void enq(int x, int y);

#define top1           q[0][head]
#define top2           q[1][head]
#define popq           tail--
#define NOMATCH       -1
#define ALL_MATCHED   -2

FILE *fp;

/******************************************************************************************/


char *getFileExtension(char *fname) {
	return strrchr(fname,'.') + 1;
}

void my_sscanf(char *in_string, char arg_decomposed[][ARG_LENGTH],int *arg_num)
{
  int i=0,j=0;
  char buffer[1024];

  *arg_num = 0;
  while(in_string[i] != 0)
  {
    if(in_string[i]!=' ' && in_string[i]!='\t' &&
      in_string[i]!='\n' && in_string[i]!=0)
    {
      buffer[j]=in_string[i];
      j++;
    }
    else if(j>0)
    {
      buffer[j] = (char)0;
      strcpy(arg_decomposed[*arg_num],buffer);
      (*arg_num)++;
      j=0;
    }
    i++;
  }
  if(j>0) /* to catch the last arg */
  {
    buffer[j] = (char)0;
    strcpy(arg_decomposed[*arg_num],buffer);
    (*arg_num)++;
  }
}

int is_number(char *in_string)
{
  int i=0;
  char c;
  while( (c=in_string[i++]) != 0)
    if((c < '0' || c >'9') && c!='.' && c!='E' && c!='e' && c!='-')
      return(0);
  return(1);
}

/***************************************************************************************/

/*********************************************************************/


int read_input(int number_case, char *fname, char *shape_file)
{
    int i,j;
    
	double x[MAXPOINTS],y[MAXPOINTS],z[MAXPOINTS];
	
    inMol = readMolecule(fname, NULL);
    mol =   babel2Mol(inMol, TRUE);
    count = mol->_size;

    if (number_case == SHAPE) {
		inShape = readMolecule(shape_file, NULL);
		shape = babel2Mol(inShape, TRUE);	
    }   

    if (permsFile != NULL) {
	FILE *permfile = fopen(permsFile,"r");
	int *used = (int *)malloc(sizeof(int) * count);		
	fscanf(permfile, "%d", &numPermutations);
	permutations = (int **)malloc(numPermutations * sizeof(int *));
	for (j = 0; j < numPermutations; j++) {		
		for (i = 0; i < count; i++) {
			used[i] = FALSE;
		}
		permutations[j] = (int *)malloc(sizeof(int) * count);	
		for (i = 0; i < count; i++) {
			int cur = -1;
			int res = fscanf(permfile,"%d", &cur);
			if (res != 1 || cur < 1 || cur > count || used[cur - 1]) {
				printf("Invalid permutation file\n");			
				free(used);
				fclose(permfile);
				exit(1);
			}
			used[cur - 1] = TRUE;
			permutations[j][i] = cur - 1;
		}
	
	}
	free(used);
	fclose(permfile);
    }    

//	printMoleculeBasic(mol);
//	printMoleculeBasic(shape);

    /*if (fscanf(fp,"%d",&count)!=1)
       { printf("EMPTY INPUT FILE\n"); exit(-1); }*/

    if(!((count == 4)||(count == 5)||(count == 6)||(count == 7)||
         (count == 8)||(count == 9)||(count ==12)||(count ==13)||
         (count ==20)||(count ==21)||(count ==60)||(count ==61) || 
	 (number_case == SHAPE && count == shape->_size)))
         { printf("ERR* NUMBER OF POINTS OUT OF RANGE %i *ERR\n", count); exit(-1); }

    switch(number_case)
       {
       case SHAPE:
       if (count != shape->_size) { printf("ERR* NUMBER VERTICES SHOULD BE: %d *ERR\n", shape->_size); }
       break;
       case B4:
       if(count!=4) { printf("ERR* NUMBER VERTICES SHOULD BE: 4 *ERR\n"); exit(1); }
       break;

       case AB4:
       if(count!=5) { printf("ERR* NUMBER VERTICES SHOULD BE: 5 *ERR\n"); exit(1); }
       break;

       case B5:
       if(count!=5) { printf("ERR* NUMBER VERTICES SHOULD BE: 5 *ERR\n"); exit(1); }
       break;

       case B5_:
       if(count!=5) { printf("ERR* NUMBER VERTICES SHOULD BE: 4 *ERR\n"); exit(1); }
       break;

       case AB5:
       if(count!=6) { printf("ERR* NUMBER VERTICES SHOULD BE: 6 *ERR\n"); exit(1); }
       break;

       case AB5_:
       if(count!=6) { printf("ERR* NUMBER VERTICES SHOULD BE: 6 *ERR\n"); exit(1); }
       break;

       case PR:
       if(count!=6) { printf("ERR* NUMBER VERTICES SHOULD BE: 6 *ERR\n"); exit(1); }
       break;

       case D3H:
       if(count!=6) { printf("ERR* NUMBER VERTICES SHOULD BE: 6 *ERR\n"); exit(1); }
       break;

       case APR:
       if(count!=7) { printf("ERR* NUMBER VERTICES SHOULD BE: 7 *ERR\n"); exit(1); }
       break;

       case AD3H:
       if(count!=7) { printf("ERR* NUMBER VERTICES SHOULD BE: 7 *ERR\n"); exit(1); }
       break;

       case PR_EQ:
       if(count!=6) { printf("ERR* NUMBER VERTICES SHOULD BE: 6 *ERR\n"); exit(1); }
       break;

       case APR_EQ:
       if(count!=7) { printf("ERR* NUMBER VERTICES SHOULD BE: 7 *ERR\n"); exit(1); }
       break;

       case B6:
       if(count!=6) { printf("ERR* NUMBER VERTICES SHOULD BE: 6 *ERR\n"); exit(1); }
       break;

       case AB6:
       if(count!=7) { printf("ERR* NUMBER VERTICES SHOULD BE: 7 *ERR\n"); exit(1); }
       break;

       case B8:
       if(count!=8) { printf("ERR* NUMBER VERTICES SHOULD BE: 8 *ERR\n"); exit(1); }
       break;

       case AB8:
       if(count!=9) { printf("ERR* NUMBER VERTICES SHOULD BE: 9 *ERR\n"); exit(1); }
       break;

       case B12:
       if(count!=12) { printf("ERR* NUMBER VERTICES SHOULD BE: 12 *ERR\n"); exit(1); }
       break;

       case AB12:
       if(count!=13) { printf("ERR* NUMBER VERTICES SHOULD BE: 13 *ERR\n"); exit(1); }
       break;

       case B20:
       if(count!=20) { printf("ERR* NUMBER VERTICES SHOULD BE: 20 *ERR\n"); exit(1); }
       break;

       case AB20:
       if(count!=21) { printf("ERR* NUMBER VERTICES SHOULD BE: 21 *ERR\n"); exit(1); }
       break;

       case B60:
       if(count!=60) { printf("ERR* NUMBER VERTICES SHOULD BE: 60 *ERR\n"); exit(1); }
       break;

       case AB60:
       if(count!=61) { printf("ERR* NUMBER VERTICES SHOULD BE: 61 *ERR\n"); exit(1); }
       break;

       case SQ5:
       if(count!=5) { printf("ERR* NUMBER VERTICES SHOULD BE: 5 *ERR\n"); exit(1); }
       break;

       default:
       printf("ERR* WRONG CASE_FLAG %i *ERR\n", number_case);
       exit(1);
       break;
       }


    for (i=0; i<count; i++)
       {
       parray[i] = alloc_matrix(3,1);
	   parray[i][1][1]=mol->_pos[i][0];
	   parray[i][2][1]=mol->_pos[i][1];
	   parray[i][3][1]=mol->_pos[i][2];

       enq(i,i);
       }

    if(number_case == 15)
       parray[0][1][1] += 0.000001;

    return(1);
}



/***********************************************************************/

void enq(int x,int y)
{
    tail++;
    q[0][tail] = x;
    q[1][tail] = y;
}

/************************************************************************/

void normalize_array(int n)
{
    double tmp;
    int i;

    x_avg = y_avg = z_avg = 0.0;
    for(i=0;i<n;i++)
       {
       x_avg += parray[i][1][1];
       y_avg += parray[i][2][1];
       z_avg += parray[i][3][1];
       }
    x_avg /= (double)(n);
    y_avg /= (double)(n);
    z_avg /= (double)(n);

    norm = 0.0;
    for(i=0;i<n;i++)
       {
       tmp = SQR(parray[i][1][1]-x_avg) +
             SQR(parray[i][2][1]-y_avg) +
             SQR(parray[i][3][1]-z_avg);
       norm += tmp;
       }
    norm = sqrt(norm/n);

    for(i=0;i<n;i++)
       {
        parray[i][1][1] = ((parray[i][1][1] - x_avg) / norm);
        parray[i][2][1] = ((parray[i][2][1] - y_avg) / norm);
        parray[i][3][1] = ((parray[i][3][1] - z_avg) / norm);
       }
}

/******************************************************************************/

void normalize_cbarray(int n) {
    double tmp, x_a, y_a, z_a, normal;
    int i;

    x_a = y_a = z_a = 0.0;
    for(i=0;i<n;i++)
       {
       x_a += cbarray[i][1][1];
       y_a += cbarray[i][2][1];
       z_a += cbarray[i][3][1];
       }
    x_a /= (double)(n);
    y_a /= (double)(n);
    z_a /= (double)(n);

    normal = 0.0;
    for(i=0;i<n;i++)
       {
       tmp = SQR(cbarray[i][1][1]-x_a) +
             SQR(cbarray[i][2][1]-y_a) +
             SQR(cbarray[i][3][1]-z_a);
       normal += tmp;
       }
    normal = sqrt(normal/n);

    for(i=0;i<n;i++)
       {
        cbarray[i][1][1] = ((cbarray[i][1][1] - x_a) / normal);
        cbarray[i][2][1] = ((cbarray[i][2][1] - y_a) / normal);
        cbarray[i][3][1] = ((cbarray[i][3][1] - z_a) / normal);
       }
}



/******************************************************************************/


void create_hedron_array(int n, int number_case) {
    int i;
    double a, b, c, d, e, f, g, h, o, p, r, t;


    for (i=0; i< n; i++)
        cbarray[i] = alloc_matrix(3,1);


/*******************************/
/*    SYMMETRIC POLYHEDRON     */
/*******************************/


    switch(number_case)
       {
       case SHAPE: {
	 int i;
	 for (i = 0; i < shape->_size; i++) {
		cbarray[i][1][1] = shape->_pos[i][0];
		cbarray[i][2][1] = shape->_pos[i][1];
		cbarray[i][3][1] = shape->_pos[i][2];		
	 }
	 break;
       }
       case B4: /*   TETRAHEDRON   */
                /*   -----------   */

          edge = sqrt(8.0/3.0);
          a    = sqrt(8.0)/3.0;
          b    = 1.0/3.0;
          c    = sqrt(2.0/3.0);
          d    = sqrt(2.0)/3.0;
          cbarray[0][1][1]=0.; cbarray[0][2][1]=0.; cbarray[0][3][1]=1.;
          cbarray[1][1][1]=0.; cbarray[1][2][1]= a; cbarray[1][3][1]=-b;
          cbarray[2][1][1]= c; cbarray[2][2][1]=-d; cbarray[2][3][1]=-b;
          cbarray[3][1][1]=-c; cbarray[3][2][1]=-d; cbarray[3][3][1]=-b;

          ADJ(0,1) = ADJ(1,0) = ADJ(0,2) = ADJ(2,0) =
          ADJ(3,0) = ADJ(0,3) = ADJ(1,2) = ADJ(2,1) =
          ADJ(1,3) = ADJ(3,1) = ADJ(2,3) = ADJ(3,2) = 1;
       break;

       case AB4: /*   TETRAHEDRON WITH CENTRE   */
                 /*   -----------------------   */

          edge = sqrt(8.0/3.0);
          a    = sqrt(8.0)/3.0;
          b    = 1.0/3.0;
          c    = sqrt(2.0/3.0);
          d    = sqrt(2.0)/3.0;
          cbarray[0][1][1]=0.; cbarray[0][2][1]=0.; cbarray[0][3][1]=1.;
          cbarray[1][1][1]=0.; cbarray[1][2][1]= a; cbarray[1][3][1]=-b;
          cbarray[2][1][1]= c; cbarray[2][2][1]=-d; cbarray[2][3][1]=-b;
          cbarray[3][1][1]=-c; cbarray[3][2][1]=-d; cbarray[3][3][1]=-b;
          cbarray[4][1][1]=0.; cbarray[4][2][1]=0.; cbarray[4][3][1]=0.;

          ADJ(0,1) = ADJ(1,0) = ADJ(0,2) = ADJ(2,0) =
          ADJ(0,3) = ADJ(3,0) = ADJ(1,2) = ADJ(2,1) =
          ADJ(1,3) = ADJ(3,1) = ADJ(2,3) = ADJ(3,2) = 1;
       break;

       case B5:  /*   BIPYRAMIDE   */
                 /*   ----------   */

          a    = sqrt(3.0/8.0);
          b    = 1.0/sqrt(8.0);
          c    = 1.0/sqrt(2.0);
          cbarray[0][1][1]=0.; cbarray[0][2][1]=0.; cbarray[0][3][1]= 1.;
          cbarray[1][1][1]=-a; cbarray[1][2][1]=-b; cbarray[1][3][1]= 0.;
          cbarray[2][1][1]= a; cbarray[2][2][1]=-b; cbarray[2][3][1]= 0.;
          cbarray[3][1][1]=0.; cbarray[3][2][1]= c; cbarray[3][3][1]= 0.;
          cbarray[4][1][1]=0.; cbarray[4][2][1]=0.; cbarray[4][3][1]=-1.;

          ADJ(0,1) = ADJ(1,0) = ADJ(0,2) = ADJ(2,0) =
          ADJ(0,3) = ADJ(3,0) = ADJ(1,2) = ADJ(2,1) =
          ADJ(1,3) = ADJ(3,1) = ADJ(1,4) = ADJ(4,1) =
          ADJ(2,3) = ADJ(3,2) = ADJ(2,4) = ADJ(4,2) =
          ADJ(3,4) = ADJ(4,3) = 1;
       break;

       case B5_:  /*   BIPYRAMIDE  (EQUIDISTANCE)   */
                  /*   --------------------------   */

          a    = sqrt(3.0)/2.0;
          b    = 0.5;
          cbarray[0][1][1]=0.; cbarray[0][2][1]=0.; cbarray[0][3][1]= 1.;
          cbarray[1][1][1]=-a; cbarray[1][2][1]=-b; cbarray[1][3][1]= 0.;
          cbarray[2][1][1]= a; cbarray[2][2][1]=-b; cbarray[2][3][1]= 0.;
          cbarray[3][1][1]=0.; cbarray[3][2][1]=1.; cbarray[3][3][1]= 0.;
          cbarray[4][1][1]=0.; cbarray[4][2][1]=0.; cbarray[4][3][1]=-1.;

          ADJ(0,1) = ADJ(1,0) = ADJ(0,2) = ADJ(2,0) =
          ADJ(0,3) = ADJ(3,0) = ADJ(1,2) = ADJ(2,1) =
          ADJ(1,3) = ADJ(3,1) = ADJ(1,4) = ADJ(4,1) =
          ADJ(2,3) = ADJ(3,2) = ADJ(2,4) = ADJ(4,2) =
          ADJ(3,4) = ADJ(4,3) = 1;
       break;


       case AB5:  /*   BIPYRAMIDE  WITH CENTRE   */
                  /*   -----------------------   */

          a    = sqrt(3.0/8.0);
          b    = 1.0/sqrt(8.0);
          c    = 1.0/sqrt(2.0);
          cbarray[0][1][1]=0.; cbarray[0][2][1]=0.; cbarray[0][3][1]= 1.;
          cbarray[1][1][1]=-a; cbarray[1][2][1]=-b; cbarray[1][3][1]= 0.;
          cbarray[2][1][1]= a; cbarray[2][2][1]=-b; cbarray[2][3][1]= 0.;
          cbarray[3][1][1]=0.; cbarray[3][2][1]= c; cbarray[3][3][1]= 0.;
          cbarray[4][1][1]=0.; cbarray[4][2][1]=0.; cbarray[4][3][1]=-1.;
          cbarray[5][1][1]=0.; cbarray[5][2][1]=0.; cbarray[5][3][1]= 0.;

          ADJ(0,1) = ADJ(1,0) = ADJ(0,2) = ADJ(2,0) =
          ADJ(0,3) = ADJ(3,0) = ADJ(1,2) = ADJ(2,1) =
          ADJ(1,3) = ADJ(3,1) = ADJ(1,4) = ADJ(4,1) =
          ADJ(2,3) = ADJ(3,2) = ADJ(2,4) = ADJ(4,2) =
          ADJ(3,4) = ADJ(4,3) = 1;
       break;


       case AB5_:  /*   BIPYRAMIDE  WITH CENTRE (EQUIDISTANCE)  */
                   /*   -------------------------------------   */

          a    = sqrt(3.0)/2.0;
          b    = 0.5;
          cbarray[0][1][1]=0.; cbarray[0][2][1]=0.; cbarray[0][3][1]= 1.;
          cbarray[1][1][1]=-a; cbarray[1][2][1]=-b; cbarray[1][3][1]= 0.;
          cbarray[2][1][1]= a; cbarray[2][2][1]=-b; cbarray[2][3][1]= 0.;
          cbarray[3][1][1]=0.; cbarray[3][2][1]=1.; cbarray[3][3][1]= 0.;
          cbarray[4][1][1]=0.; cbarray[4][2][1]=0.; cbarray[4][3][1]=-1.;
          cbarray[5][1][1]=0.; cbarray[5][2][1]=0.; cbarray[5][3][1]= 0.;

          ADJ(0,1) = ADJ(1,0) = ADJ(0,2) = ADJ(2,0) =
          ADJ(0,3) = ADJ(3,0) = ADJ(1,2) = ADJ(2,1) =
          ADJ(1,3) = ADJ(3,1) = ADJ(1,4) = ADJ(4,1) =
          ADJ(2,3) = ADJ(3,2) = ADJ(2,4) = ADJ(4,2) =
          ADJ(3,4) = ADJ(4,3) = 1;
       break;

       case PR:  /*   TRIGONAL PRISM   */
                 /*   --------------   */

          a = 1.0/sqrt(2.0);
          b = 1.0/sqrt(3.0);
          c = 1.0/sqrt(6.0);
          d = 2.0*c;

          cbarray[0][1][1]=0.; cbarray[0][2][1]=-d; cbarray[0][3][1]= b;
          cbarray[1][1][1]=-a; cbarray[1][2][1]= c; cbarray[1][3][1]= b;
          cbarray[2][1][1]= a; cbarray[2][2][1]= c; cbarray[2][3][1]= b;
          cbarray[3][1][1]=0.; cbarray[3][2][1]=-d; cbarray[3][3][1]=-b;
          cbarray[4][1][1]=-a; cbarray[4][2][1]= c; cbarray[4][3][1]=-b;
          cbarray[5][1][1]= a; cbarray[5][2][1]= c; cbarray[5][3][1]=-b;

          ADJ(0,1) = ADJ(1,0) = ADJ(0,2) = ADJ(2,0) =
          ADJ(0,3) = ADJ(3,0) = ADJ(1,2) = ADJ(2,1) =
          ADJ(1,4) = ADJ(4,1) = ADJ(2,5) = ADJ(5,2) =
          ADJ(3,4) = ADJ(4,3) = ADJ(3,5) = ADJ(5,3) =
          ADJ(4,5) = ADJ(5,4) = 1;
       break;

       case APR:  /*   TRIGONAL PRISM WITH CENTRE  */
                  /*   --------------------------  */

          a = 1.0/sqrt(2.0);
          b = 1.0/sqrt(3.0);
          c = 1.0/sqrt(6.0);
          d = 2.0*c;

          cbarray[0][1][1]=0.; cbarray[0][2][1]=-d; cbarray[0][3][1]= b;
          cbarray[1][1][1]=-a; cbarray[1][2][1]= c; cbarray[1][3][1]= b;
          cbarray[2][1][1]= a; cbarray[2][2][1]= c; cbarray[2][3][1]= b;
          cbarray[3][1][1]=0.; cbarray[3][2][1]=-d; cbarray[3][3][1]=-b;
          cbarray[4][1][1]=-a; cbarray[4][2][1]= c; cbarray[4][3][1]=-b;
          cbarray[5][1][1]= a; cbarray[5][2][1]= c; cbarray[5][3][1]=-b;
          cbarray[6][1][1]=0.; cbarray[6][2][1]=0.; cbarray[6][3][1]=0.;

          ADJ(0,1) = ADJ(1,0) = ADJ(0,2) = ADJ(2,0) =
          ADJ(0,3) = ADJ(3,0) = ADJ(1,2) = ADJ(2,1) =
          ADJ(1,4) = ADJ(4,1) = ADJ(2,5) = ADJ(5,2) =
          ADJ(3,4) = ADJ(4,3) = ADJ(3,5) = ADJ(5,3) =
          ADJ(4,5) = ADJ(5,4) = 1;
       break;

       case D3H:  /*   TRIGONAL PRISM (D3h - symmetry)   */
                  /*   -------------------------------   */

          a = 1.0/sqrt(2.0);
          b = 1.0/sqrt(3.0);
          c = 1.0/sqrt(6.0);
          d = 2.0*c;

          cbarray[0][1][1]=0.; cbarray[0][2][1]=-d; cbarray[0][3][1]= b;
          cbarray[1][1][1]=-a; cbarray[1][2][1]= c; cbarray[1][3][1]= b;
          cbarray[2][1][1]= a; cbarray[2][2][1]= c; cbarray[2][3][1]= b;
          cbarray[3][1][1]=0.; cbarray[3][2][1]=-d; cbarray[3][3][1]=-b;
          cbarray[4][1][1]=-a; cbarray[4][2][1]= c; cbarray[4][3][1]=-b;
          cbarray[5][1][1]= a; cbarray[5][2][1]= c; cbarray[5][3][1]=-b;

          ADJ(0,1) = ADJ(1,0) = ADJ(0,2) = ADJ(2,0) =
          ADJ(0,3) = ADJ(3,0) = ADJ(1,2) = ADJ(2,1) =
          ADJ(1,4) = ADJ(4,1) = ADJ(2,5) = ADJ(5,2) =
          ADJ(3,4) = ADJ(4,3) = ADJ(3,5) = ADJ(5,3) =
          ADJ(4,5) = ADJ(5,4) = 1;
       break;

       case AD3H:  /*   TRIGONAL PRISM WITH CENTRE (D3h - symmetry) */
                   /*   ------------------------------------------  */

          a = 1.0/sqrt(2.0);
          b = 1.0/sqrt(3.0);
          c = 1.0/sqrt(6.0);
          d = 2.0*c;

          cbarray[0][1][1]=0.; cbarray[0][2][1]=-d; cbarray[0][3][1]= b;
          cbarray[1][1][1]=-a; cbarray[1][2][1]= c; cbarray[1][3][1]= b;
          cbarray[2][1][1]= a; cbarray[2][2][1]= c; cbarray[2][3][1]= b;
          cbarray[3][1][1]=0.; cbarray[3][2][1]=-d; cbarray[3][3][1]=-b;
          cbarray[4][1][1]=-a; cbarray[4][2][1]= c; cbarray[4][3][1]=-b;
          cbarray[5][1][1]= a; cbarray[5][2][1]= c; cbarray[5][3][1]=-b;
          cbarray[6][1][1]=0.; cbarray[6][2][1]=0.; cbarray[6][3][1]=0.;

          ADJ(0,1) = ADJ(1,0) = ADJ(0,2) = ADJ(2,0) =
          ADJ(0,3) = ADJ(3,0) = ADJ(1,2) = ADJ(2,1) =
          ADJ(1,4) = ADJ(4,1) = ADJ(2,5) = ADJ(5,2) =
          ADJ(3,4) = ADJ(4,3) = ADJ(3,5) = ADJ(5,3) =
          ADJ(4,5) = ADJ(5,4) = 1;
       break;

       case PR_EQ:  /*   TRIGONAL EQUILATERAL PRISM   */
                    /*   --------------------------   */

          a = 1.0/sqrt(2.0);
          b = 1.0/sqrt(3.0);
          c = 1.0/sqrt(6.0);
          d = 2.0*c;

          cbarray[0][1][1]=0.; cbarray[0][2][1]=-d; cbarray[0][3][1]= a;
          cbarray[1][1][1]=-a; cbarray[1][2][1]= c; cbarray[1][3][1]= a;
          cbarray[2][1][1]= a; cbarray[2][2][1]= c; cbarray[2][3][1]= a;
          cbarray[3][1][1]=0.; cbarray[3][2][1]=-d; cbarray[3][3][1]=-a;
          cbarray[4][1][1]=-a; cbarray[4][2][1]= c; cbarray[4][3][1]=-a;
          cbarray[5][1][1]= a; cbarray[5][2][1]= c; cbarray[5][3][1]=-a;

          ADJ(0,1) = ADJ(1,0) = ADJ(0,2) = ADJ(2,0) =
          ADJ(0,3) = ADJ(3,0) = ADJ(1,2) = ADJ(2,1) =
          ADJ(1,4) = ADJ(4,1) = ADJ(2,5) = ADJ(5,2) =
          ADJ(3,4) = ADJ(4,3) = ADJ(3,5) = ADJ(5,3) =
          ADJ(4,5) = ADJ(5,4) = 1;
       break;

       case APR_EQ:  /*   TRIGONAL EQUILATERAL PRISM WITH CENTRE  */
                     /*   --------------------------------------  */

          a = 1.0/sqrt(2.0);
          b = 1.0/sqrt(3.0);
          c = 1.0/sqrt(6.0);
          d = 2.0*c;

          cbarray[0][1][1]=0.; cbarray[0][2][1]=-d; cbarray[0][3][1]= a;
          cbarray[1][1][1]=-a; cbarray[1][2][1]= c; cbarray[1][3][1]= a;
          cbarray[2][1][1]= a; cbarray[2][2][1]= c; cbarray[2][3][1]= a;
          cbarray[3][1][1]=0.; cbarray[3][2][1]=-d; cbarray[3][3][1]=-a;
          cbarray[4][1][1]=-a; cbarray[4][2][1]= c; cbarray[4][3][1]=-a;
          cbarray[5][1][1]= a; cbarray[5][2][1]= c; cbarray[5][3][1]=-a;
          cbarray[6][1][1]=0.; cbarray[6][2][1]=0.; cbarray[6][3][1]=0.;

          ADJ(0,1) = ADJ(1,0) = ADJ(0,2) = ADJ(2,0) =
          ADJ(0,3) = ADJ(3,0) = ADJ(1,2) = ADJ(2,1) =
          ADJ(1,4) = ADJ(4,1) = ADJ(2,5) = ADJ(5,2) =
          ADJ(3,4) = ADJ(4,3) = ADJ(3,5) = ADJ(5,3) =
          ADJ(4,5) = ADJ(5,4) = 1;
       break;



       case B6:  /*   OCTAHEDRON   */
                 /*   ----------   */

          edge = sqrt(2.0);
          a    = 1.0/sqrt(2.0);
          cbarray[0][1][1]=0.; cbarray[0][2][1]=0.; cbarray[0][3][1]= 1.;
          cbarray[1][1][1]= a; cbarray[1][2][1]= a; cbarray[1][3][1]= 0.;
          cbarray[2][1][1]=-a; cbarray[2][2][1]= a; cbarray[2][3][1]= 0.;
          cbarray[3][1][1]=-a; cbarray[3][2][1]=-a; cbarray[3][3][1]= 0.;
          cbarray[4][1][1]= a; cbarray[4][2][1]=-a; cbarray[4][3][1]= 0.;
          cbarray[5][1][1]=0.; cbarray[5][2][1]=0.; cbarray[5][3][1]=-1.;

          ADJ(0,1) = ADJ(1,0) = ADJ(0,2) = ADJ(2,0) =
          ADJ(0,3) = ADJ(3,0) = ADJ(0,4) = ADJ(4,0) =
          ADJ(1,2) = ADJ(2,1) = ADJ(1,4) = ADJ(4,1) =
          ADJ(1,5) = ADJ(5,1) = ADJ(2,3) = ADJ(3,2) =
          ADJ(2,5) = ADJ(5,2) = ADJ(3,4) = ADJ(4,3) =
          ADJ(3,5) = ADJ(5,3) = ADJ(4,5) = ADJ(5,4) = 1;
       break;

       case AB6:  /*   OCTAHEDRON WITH CENTRE   */
                  /*   ----------------------   */
          edge = sqrt(2.0);
          a    = 1.0/sqrt(2.0);
          cbarray[0][1][1]=0.; cbarray[0][2][1]=0.; cbarray[0][3][1]= 1.;
          cbarray[1][1][1]= a; cbarray[1][2][1]= a; cbarray[1][3][1]= 0.;
          cbarray[2][1][1]=-a; cbarray[2][2][1]= a; cbarray[2][3][1]= 0.;
          cbarray[3][1][1]=-a; cbarray[3][2][1]=-a; cbarray[3][3][1]= 0.;
          cbarray[4][1][1]= a; cbarray[4][2][1]=-a; cbarray[4][3][1]= 0.;
          cbarray[5][1][1]=0.; cbarray[5][2][1]=0.; cbarray[5][3][1]=-1.;
          cbarray[6][1][1]=0.; cbarray[6][2][1]=0.; cbarray[6][3][1]= 0.;

          ADJ(0,1) = ADJ(1,0) = ADJ(0,2) = ADJ(2,0) =
          ADJ(0,3) = ADJ(3,0) = ADJ(0,4) = ADJ(4,0) =
          ADJ(1,2) = ADJ(2,1) = ADJ(1,4) = ADJ(4,1) =
          ADJ(1,5) = ADJ(5,1) = ADJ(2,3) = ADJ(3,2) =
          ADJ(2,5) = ADJ(5,2) = ADJ(3,4) = ADJ(4,3) =
          ADJ(3,5) = ADJ(5,3) = ADJ(4,5) = ADJ(5,4) = 1;
       break;

       case B8:  /*   CUBE   */
                 /*   ----   */

          edge = 2.0/sqrt(3.0);
          a    = 1.0/sqrt(3.0);
          cbarray[0][1][1]= a; cbarray[0][2][1]= a; cbarray[0][3][1]= a;
          cbarray[1][1][1]=-a; cbarray[1][2][1]= a; cbarray[1][3][1]= a;
          cbarray[2][1][1]= a; cbarray[2][2][1]=-a; cbarray[2][3][1]= a;
          cbarray[3][1][1]= a; cbarray[3][2][1]= a; cbarray[3][3][1]=-a;
          cbarray[4][1][1]=-a; cbarray[4][2][1]=-a; cbarray[4][3][1]= a;
          cbarray[5][1][1]= a; cbarray[5][2][1]=-a; cbarray[5][3][1]=-a;
          cbarray[6][1][1]=-a; cbarray[6][2][1]= a; cbarray[6][3][1]=-a;
          cbarray[7][1][1]=-a; cbarray[7][2][1]=-a; cbarray[7][3][1]=-a;

          ADJ(0,1) = ADJ(1,0) = ADJ(0,2) = ADJ(2,0) =
          ADJ(0,3) = ADJ(3,0) = ADJ(1,4) = ADJ(4,1) =
          ADJ(1,6) = ADJ(6,1) = ADJ(2,4) = ADJ(4,2) =
          ADJ(2,5) = ADJ(5,2) = ADJ(3,5) = ADJ(5,3) =
          ADJ(3,6) = ADJ(6,3) = ADJ(4,7) = ADJ(7,4) =
          ADJ(5,7) = ADJ(7,5) = ADJ(6,7) = ADJ(7,6) = 1;
       break;

       case AB8:  /*   CUBE WITH CENTRE   */
                  /*   ----------------   */

          edge = 2.0/sqrt(3.0);
          a    = 1.0/sqrt(3.0);
          cbarray[0][1][1]= a; cbarray[0][2][1]= a; cbarray[0][3][1]= a;
          cbarray[1][1][1]=-a; cbarray[1][2][1]= a; cbarray[1][3][1]= a;
          cbarray[2][1][1]= a; cbarray[2][2][1]=-a; cbarray[2][3][1]= a;
          cbarray[3][1][1]= a; cbarray[3][2][1]= a; cbarray[3][3][1]=-a;
          cbarray[4][1][1]=-a; cbarray[4][2][1]=-a; cbarray[4][3][1]= a;
          cbarray[5][1][1]= a; cbarray[5][2][1]=-a; cbarray[5][3][1]=-a;
          cbarray[6][1][1]=-a; cbarray[6][2][1]= a; cbarray[6][3][1]=-a;
          cbarray[7][1][1]=-a; cbarray[7][2][1]=-a; cbarray[7][3][1]=-a;
          cbarray[8][1][1]=0.; cbarray[8][2][1]=0.; cbarray[8][3][1]=0.;

          ADJ(0,1) = ADJ(1,0) = ADJ(0,2) = ADJ(2,0) =
          ADJ(0,3) = ADJ(3,0) = ADJ(1,4) = ADJ(4,1) =
          ADJ(1,6) = ADJ(6,1) = ADJ(2,4) = ADJ(4,2) =
          ADJ(2,5) = ADJ(5,2) = ADJ(3,5) = ADJ(5,3) =
          ADJ(3,6) = ADJ(6,3) = ADJ(4,7) = ADJ(7,4) =
          ADJ(5,7) = ADJ(7,5) = ADJ(6,7) = ADJ(7,6) = 1;
       break;

       case B12:  /*   ICOSAHEDRON   */
                  /*   -----------   */

          edge = sqrt(8.0/(5.0+sqrt(5.0)));
          r    = edge*sqrt(1.0-edge*edge/4.0);
          h    = 1.0-edge*edge/2.0;
          a    = r*sin(TWOPI/5.0);
          b    = r*cos(TWOPI/5.0);
          c    = r*sin(PI/5.0);
          d    = r*cos(PI/5.0);
          cbarray[ 0][1][1]=0.; cbarray[ 0][2][1]=0.; cbarray[ 0][3][1]=1.;
          cbarray[ 1][1][1]=0.; cbarray[ 1][2][1]=-r; cbarray[ 1][3][1]= h;
          cbarray[ 2][1][1]=-a; cbarray[ 2][2][1]=-b; cbarray[ 2][3][1]= h;
          cbarray[ 3][1][1]=-c; cbarray[ 3][2][1]= d; cbarray[ 3][3][1]= h;
          cbarray[ 4][1][1]= c; cbarray[ 4][2][1]= d; cbarray[ 4][3][1]= h;
          cbarray[ 5][1][1]= a; cbarray[ 5][2][1]=-b; cbarray[ 5][3][1]= h;
          cbarray[ 6][1][1]=-c; cbarray[ 6][2][1]=-d; cbarray[ 6][3][1]=-h;
          cbarray[ 7][1][1]=-a; cbarray[ 7][2][1]= b; cbarray[ 7][3][1]=-h;
          cbarray[ 8][1][1]=0.; cbarray[ 8][2][1]= r; cbarray[ 8][3][1]=-h;
          cbarray[ 9][1][1]= a; cbarray[ 9][2][1]= b; cbarray[ 9][3][1]=-h;
          cbarray[10][1][1]= c; cbarray[10][2][1]=-d; cbarray[10][3][1]=-h;
          cbarray[11][1][1]=0.; cbarray[11][2][1]=0.; cbarray[11][3][1]=-1.;

          ADJ(0, 1) = ADJ(1 ,0) = ADJ(0 , 2) = ADJ(2 , 0) =
          ADJ(0, 3) = ADJ(3 ,0) = ADJ(0 , 4) = ADJ(4 , 0) =
          ADJ(0, 5) = ADJ(5 ,0) = ADJ(1 , 2) = ADJ(2 , 1) =
          ADJ(1, 5) = ADJ(5 ,1) = ADJ(1 , 6) = ADJ(6 , 1) =
          ADJ(1,10) = ADJ(10,1) = ADJ(2 , 3) = ADJ(3 , 2) =
          ADJ(2, 6) = ADJ(6 ,2) = ADJ(2 , 7) = ADJ(7 , 2) =
          ADJ(3, 4) = ADJ(4 ,3) = ADJ(3 , 7) = ADJ(7 , 3) =
          ADJ(3, 8) = ADJ(8 ,3) = ADJ(4 , 5) = ADJ(5 , 4) =
          ADJ(4, 8) = ADJ(8 ,4) = ADJ(4 , 9) = ADJ(9 , 4) =
          ADJ(5, 9) = ADJ(9 ,5) = ADJ(5 ,10) = ADJ(10, 5) =
          ADJ(6, 7) = ADJ(7 ,6) = ADJ(6 ,10) = ADJ(10, 6) =
          ADJ(6,11) = ADJ(11,6) = ADJ(7 , 8) = ADJ(8 , 7) =
          ADJ(7,11) = ADJ(11,7) = ADJ(8 , 9) = ADJ(9 , 8) =
          ADJ(8,11) = ADJ(11,8) = ADJ(9 ,10) = ADJ(10, 9) =
          ADJ(9,11) = ADJ(11,9) = ADJ(10,11) = ADJ(11,10) = 1;
       break;

       case AB12:  /*   ICOSAHEDRON WITH CENTRE  */
                   /*   -----------------------  */

          edge = sqrt(8.0/(5.0+sqrt(5.0)));
          r    = edge*sqrt(1.0-edge*edge/4.0);
          h    = 1.0-edge*edge/2.0;
          a    = r*sin(TWOPI/5.0);
          b    = r*cos(TWOPI/5.0);
          c    = r*sin(PI/5.0);
          d    = r*cos(PI/5.0);
          cbarray[ 0][1][1]=0.; cbarray[ 0][2][1]=0.; cbarray[ 0][3][1]=1.;
          cbarray[ 1][1][1]=0.; cbarray[ 1][2][1]=-r; cbarray[ 1][3][1]= h;
          cbarray[ 2][1][1]=-a; cbarray[ 2][2][1]=-b; cbarray[ 2][3][1]= h;
          cbarray[ 3][1][1]=-c; cbarray[ 3][2][1]= d; cbarray[ 3][3][1]= h;
          cbarray[ 4][1][1]= c; cbarray[ 4][2][1]= d; cbarray[ 4][3][1]= h;
          cbarray[ 5][1][1]= a; cbarray[ 5][2][1]=-b; cbarray[ 5][3][1]= h;
          cbarray[ 6][1][1]=-c; cbarray[ 6][2][1]=-d; cbarray[ 6][3][1]=-h;
          cbarray[ 7][1][1]=-a; cbarray[ 7][2][1]= b; cbarray[ 7][3][1]=-h;
          cbarray[ 8][1][1]=0.; cbarray[ 8][2][1]= r; cbarray[ 8][3][1]=-h;
          cbarray[ 9][1][1]= a; cbarray[ 9][2][1]= b; cbarray[ 9][3][1]=-h;
          cbarray[10][1][1]= c; cbarray[10][2][1]=-d; cbarray[10][3][1]=-h;
          cbarray[11][1][1]=0.; cbarray[11][2][1]=0.; cbarray[11][3][1]=-1.;
          cbarray[12][1][1]=0.; cbarray[12][2][1]=0.; cbarray[12][3][1]=0.;

          ADJ(0, 1) = ADJ(1 ,0) = ADJ(0 , 2) = ADJ(2 , 0) =
          ADJ(0, 3) = ADJ(3 ,0) = ADJ(0 , 4) = ADJ(4 , 0) =
          ADJ(0, 5) = ADJ(5 ,0) = ADJ(1 , 2) = ADJ(2 , 1) =
          ADJ(1, 5) = ADJ(5 ,1) = ADJ(1 , 6) = ADJ(6 , 1) =
          ADJ(1,10) = ADJ(10,1) = ADJ(2 , 3) = ADJ(3 , 2) =
          ADJ(2, 6) = ADJ(6 ,2) = ADJ(2 , 7) = ADJ(7 , 2) =
          ADJ(3, 4) = ADJ(4 ,3) = ADJ(3 , 7) = ADJ(7 , 3) =
          ADJ(3, 8) = ADJ(8 ,3) = ADJ(4 , 5) = ADJ(5 , 4) =
          ADJ(4, 8) = ADJ(8 ,4) = ADJ(4 , 9) = ADJ(9 , 4) =
          ADJ(5, 9) = ADJ(9 ,5) = ADJ(5 ,10) = ADJ(10, 5) =
          ADJ(6, 7) = ADJ(7 ,6) = ADJ(6 ,10) = ADJ(10, 6) =
          ADJ(6,11) = ADJ(11,6) = ADJ(7 , 8) = ADJ(8 , 7) =
          ADJ(7,11) = ADJ(11,7) = ADJ(8 , 9) = ADJ(9 , 8) =
          ADJ(8,11) = ADJ(11,8) = ADJ(9 ,10) = ADJ(10, 9) =
          ADJ(9,11) = ADJ(11,9) = ADJ(10,11) = ADJ(11,10) = 1;
       break;

       case B20:  /*   DODECAHEDRON   */
                  /*   ------------   */

          edge = 4.0/(sqrt(3.0)*(1.0+sqrt(5.0)));
          a    = edge/(2.0*sin(PI/5.0));
          c    = sqrt(10.0+22.0/sqrt(5.0))/(sqrt(3.0)*(1.0+sqrt(5.0)));
          d    = c-(c*edge*edge+edge*a*sqrt(4.0-edge*edge))/2.0;
          b    = sqrt(1.0-d*d);
          e    = a*sin(TWOPI/5.0);
          f    = a*cos(TWOPI/5.0);
          g    = a*sin(PI/5.0);
          h    = a*cos(PI/5.0);
          o    = b*sin(TWOPI/5.0);
          p    = b*cos(TWOPI/5.0);
          r    = b*sin(PI/5.0);
          t    = b*cos(PI/5.0);
          cbarray[ 0][1][1]=0.; cbarray[ 0][2][1]=-a; cbarray[ 0][3][1]= c;
          cbarray[ 1][1][1]=-e; cbarray[ 1][2][1]=-f; cbarray[ 1][3][1]= c;
          cbarray[ 2][1][1]=-g; cbarray[ 2][2][1]= h; cbarray[ 2][3][1]= c;
          cbarray[ 3][1][1]= g; cbarray[ 3][2][1]= h; cbarray[ 3][3][1]= c;
          cbarray[ 4][1][1]= e; cbarray[ 4][2][1]=-f; cbarray[ 4][3][1]= c;
          cbarray[ 5][1][1]=0.; cbarray[ 5][2][1]=-b; cbarray[ 5][3][1]= d;
          cbarray[ 6][1][1]=-o; cbarray[ 6][2][1]=-p; cbarray[ 6][3][1]= d;
          cbarray[ 7][1][1]=-r; cbarray[ 7][2][1]= t; cbarray[ 7][3][1]= d;
          cbarray[ 8][1][1]= r; cbarray[ 8][2][1]= t; cbarray[ 8][3][1]= d;
          cbarray[ 9][1][1]= o; cbarray[ 9][2][1]=-p; cbarray[ 9][3][1]= d;
          cbarray[10][1][1]=-r; cbarray[10][2][1]=-t; cbarray[10][3][1]=-d;
          cbarray[11][1][1]= r; cbarray[11][2][1]=-t; cbarray[11][3][1]=-d;
          cbarray[12][1][1]=-o; cbarray[12][2][1]= p; cbarray[12][3][1]=-d;
          cbarray[13][1][1]=0.; cbarray[13][2][1]= b; cbarray[13][3][1]=-d;
          cbarray[14][1][1]= o; cbarray[14][2][1]= p; cbarray[14][3][1]=-d;
          cbarray[15][1][1]=-g; cbarray[15][2][1]=-h; cbarray[15][3][1]=-c;
          cbarray[16][1][1]= g; cbarray[16][2][1]=-h; cbarray[16][3][1]=-c;
          cbarray[17][1][1]=-e; cbarray[17][2][1]= f; cbarray[17][3][1]=-c;
          cbarray[18][1][1]=0.; cbarray[18][2][1]= a; cbarray[18][3][1]=-c;
          cbarray[19][1][1]= e; cbarray[19][2][1]= f; cbarray[19][3][1]=-c;

          ADJ(0 , 1) = ADJ(1 , 0) = ADJ(0 , 4) = ADJ(4 , 0) =
          ADJ(0 , 5) = ADJ(5 , 0) = ADJ(1 , 2) = ADJ(2 , 1) =
          ADJ(1 , 6) = ADJ(6 , 1) = ADJ(2 , 3) = ADJ(3 , 2) =
          ADJ(2 , 7) = ADJ(7 , 2) = ADJ(3 , 4) = ADJ(4 , 3) =
          ADJ(3 , 8) = ADJ(8 , 3) = ADJ(4 , 9) = ADJ(9 , 4) =
          ADJ(5 ,10) = ADJ(10, 5) = ADJ(5 ,11) = ADJ(11, 5) =
          ADJ(6 ,10) = ADJ(10 ,6) = ADJ(6 ,12) = ADJ(12, 6) =
          ADJ(7 ,12) = ADJ(12, 7) = ADJ(7 ,13) = ADJ(13, 7) =
          ADJ(8 ,13) = ADJ(13, 8) = ADJ(8 ,14) = ADJ(14, 8) =
          ADJ(9 ,11) = ADJ(11, 9) = ADJ(9 ,14) = ADJ(14, 9) =
          ADJ(10,15) = ADJ(15,10) = ADJ(11,16) = ADJ(16,11) =
          ADJ(12,17) = ADJ(17,12) = ADJ(13,18) = ADJ(18,13) =
          ADJ(14,19) = ADJ(19,14) = ADJ(15,16) = ADJ(16,15) =
          ADJ(15,17) = ADJ(17,15) = ADJ(16,19) = ADJ(19,16) =
          ADJ(17,18) = ADJ(18,17) = ADJ(18,19) = ADJ(19,18) = 1;
       break;

       case AB20:  /*   DODECAHEDRON WITH CENTRE   */
                   /*   ------------------------   */

          edge = 4.0/(sqrt(3.0)*(1.0+sqrt(5.0)));
          a    = edge/(2.0*sin(PI/5.0));
          c    = sqrt(10.0+22.0/sqrt(5.0))/(sqrt(3.0)*(1.0+sqrt(5.0)));
          d    = c-(c*edge*edge+edge*a*sqrt(4.0-edge*edge))/2.0;
          b    = sqrt(1.0-d*d);
          e    = a*sin(TWOPI/5.0);
          f    = a*cos(TWOPI/5.0);
          g    = a*sin(PI/5.0);
          h    = a*cos(PI/5.0);
          o    = b*sin(TWOPI/5.0);
          p    = b*cos(TWOPI/5.0);
          r    = b*sin(PI/5.0);
          t    = b*cos(PI/5.0);
          cbarray[ 0][1][1]=0.; cbarray[ 0][2][1]=-a; cbarray[ 0][3][1]= c;
          cbarray[ 1][1][1]=-e; cbarray[ 1][2][1]=-f; cbarray[ 1][3][1]= c;
          cbarray[ 2][1][1]=-g; cbarray[ 2][2][1]= h; cbarray[ 2][3][1]= c;
          cbarray[ 3][1][1]= g; cbarray[ 3][2][1]= h; cbarray[ 3][3][1]= c;
          cbarray[ 4][1][1]= e; cbarray[ 4][2][1]=-f; cbarray[ 4][3][1]= c;
          cbarray[ 5][1][1]=0.; cbarray[ 5][2][1]=-b; cbarray[ 5][3][1]= d;
          cbarray[ 6][1][1]=-o; cbarray[ 6][2][1]=-p; cbarray[ 6][3][1]= d;
          cbarray[ 7][1][1]=-r; cbarray[ 7][2][1]= t; cbarray[ 7][3][1]= d;
          cbarray[ 8][1][1]= r; cbarray[ 8][2][1]= t; cbarray[ 8][3][1]= d;
          cbarray[ 9][1][1]= o; cbarray[ 9][2][1]=-p; cbarray[ 9][3][1]= d;
          cbarray[10][1][1]=-r; cbarray[10][2][1]=-t; cbarray[10][3][1]=-d;
          cbarray[11][1][1]= r; cbarray[11][2][1]=-t; cbarray[11][3][1]=-d;
          cbarray[12][1][1]=-o; cbarray[12][2][1]= p; cbarray[12][3][1]=-d;
          cbarray[13][1][1]=0.; cbarray[13][2][1]= b; cbarray[13][3][1]=-d;
          cbarray[14][1][1]= o; cbarray[14][2][1]= p; cbarray[14][3][1]=-d;
          cbarray[15][1][1]=-g; cbarray[15][2][1]=-h; cbarray[15][3][1]=-c;
          cbarray[16][1][1]= g; cbarray[16][2][1]=-h; cbarray[16][3][1]=-c;
          cbarray[17][1][1]=-e; cbarray[17][2][1]= f; cbarray[17][3][1]=-c;
          cbarray[18][1][1]=0.; cbarray[18][2][1]= a; cbarray[18][3][1]=-c;
          cbarray[19][1][1]= e; cbarray[19][2][1]= f; cbarray[19][3][1]=-c;
          cbarray[20][1][1]=0.; cbarray[20][2][1]=0.; cbarray[20][3][1]=0.;

          ADJ(0 , 1) = ADJ(1 , 0) = ADJ(0 , 4) = ADJ(4 , 0) =
          ADJ(0 , 5) = ADJ(5 , 0) = ADJ(1 , 2) = ADJ(2 , 1) =
          ADJ(1 , 6) = ADJ(6 , 1) = ADJ(2 , 3) = ADJ(3 , 2) =
          ADJ(2 , 7) = ADJ(7 , 2) = ADJ(3 , 4) = ADJ(4 , 3) =
          ADJ(3 , 8) = ADJ(8 , 3) = ADJ(4 , 9) = ADJ(9 , 4) =
          ADJ(5 ,10) = ADJ(10, 5) = ADJ(5 ,11) = ADJ(11, 5) =
          ADJ(6 ,10) = ADJ(10 ,6) = ADJ(6 ,12) = ADJ(12, 6) =
          ADJ(7 ,12) = ADJ(12, 7) = ADJ(7 ,13) = ADJ(13, 7) =
          ADJ(8 ,13) = ADJ(13, 8) = ADJ(8 ,14) = ADJ(14, 8) =
          ADJ(9 ,11) = ADJ(11, 9) = ADJ(9 ,14) = ADJ(14, 9) =
          ADJ(10,15) = ADJ(15,10) = ADJ(11,16) = ADJ(16,11) =
          ADJ(12,17) = ADJ(17,12) = ADJ(13,18) = ADJ(18,13) =
          ADJ(14,19) = ADJ(19,14) = ADJ(15,16) = ADJ(16,15) =
          ADJ(15,17) = ADJ(17,15) = ADJ(16,19) = ADJ(19,16) =
          ADJ(17,18) = ADJ(18,17) = ADJ(18,19) = ADJ(19,18) = 1;
       break;


       case B60:  /*   THE BALL   */
                  /*   --------   */

/*
cbarray[ 0][1][1]= 0.000000; cbarray[ 0][2][1]= 0.000000; cbarray[ 0][3][1]= 6.798521;
cbarray[ 1][1][1]= 2.020076; cbarray[ 1][2][1]= 1.578943; cbarray[ 1][3][1]= 6.296516;
cbarray[ 2][1][1]= 1.069376; cbarray[ 2][2][1]= 3.872887; cbarray[ 2][3][1]= 5.484256;
cbarray[ 3][1][1]=-2.199175; cbarray[ 3][2][1]= 1.318104; cbarray[ 3][3][1]= 6.296516;
cbarray[ 4][1][1]=-1.538264; cbarray[ 4][2][1]= 3.711680; cbarray[ 4][3][1]= 5.484256;
cbarray[ 5][1][1]=-4.523978; cbarray[ 5][2][1]=-0.078197; cbarray[ 5][3][1]= 5.074189;
cbarray[ 6][1][1]=-4.366675; cbarray[ 6][2][1]=-2.622671; cbarray[ 6][3][1]= 4.502625;
cbarray[ 7][1][1]= 0.179099; cbarray[ 7][2][1]=-2.897047; cbarray[ 7][3][1]= 6.147758;
cbarray[ 8][1][1]=-1.862773; cbarray[ 8][2][1]=-4.123417; cbarray[ 8][3][1]= 5.074189;
cbarray[ 9][1][1]= 2.356477; cbarray[ 9][2][1]=-3.862577; cbarray[ 9][3][1]= 5.074189;
cbarray[10][1][1]= 4.499161; cbarray[10][2][1]= 0.479624; cbarray[10][3][1]= 5.074189;
cbarray[11][1][1]= 4.656463; cbarray[11][2][1]=-2.064849; cbarray[11][3][1]= 4.502625;
cbarray[12][1][1]= 5.725840; cbarray[12][2][1]= 1.808038; cbarray[12][3][1]= 3.188360;
cbarray[13][1][1]= 2.466029; cbarray[13][2][1]= 5.385371; cbarray[13][3][1]= 3.337118;
cbarray[14][1][1]= 4.643408; cbarray[14][2][1]= 4.419841; cbarray[14][3][1]= 2.263549;
cbarray[15][1][1]= 1.085068; cbarray[15][2][1]= 6.552577; cbarray[15][3][1]= 1.451289;
cbarray[16][1][1]=-3.110577; cbarray[16][2][1]= 5.040619; cbarray[16][3][1]= 3.337118;
cbarray[17][1][1]=-1.883898; cbarray[17][2][1]= 6.369032; cbarray[17][3][1]= 1.451289;
cbarray[18][1][1]=-5.904939; cbarray[18][2][1]= 1.089009; cbarray[18][3][1]= 3.188360;
cbarray[19][1][1]=-5.152449; cbarray[19][2][1]= 3.814249; cbarray[19][3][1]= 2.263549;
cbarray[20][1][1]=-5.650418; cbarray[20][2][1]=-3.028036; cbarray[20][3][1]= 2.263549;
cbarray[21][1][1]=-6.601117; cbarray[21][2][1]=-0.734092; cbarray[21][3][1]= 1.451289;
cbarray[22][1][1]=-4.608140; cbarray[22][2][1]=-4.990316; cbarray[22][3][1]= 0.285782;
cbarray[23][1][1]=-0.947341; cbarray[23][2][1]=-5.846885; cbarray[23][3][1]= 3.337118;
cbarray[24][1][1]= 1.660299; cbarray[24][2][1]=-5.685678; cbarray[24][3][1]= 3.337118;
cbarray[25][1][1]=-2.408965; cbarray[25][2][1]=-6.308420; cbarray[25][3][1]= 0.787787;
cbarray[26][1][1]= 3.167641; cbarray[26][2][1]=-5.963667; cbarray[26][3][1]= 0.787787;
cbarray[27][1][1]= 5.980361; cbarray[27][2][1]=-2.309007; cbarray[27][3][1]= 2.263549;
cbarray[28][1][1]= 6.641272; cbarray[28][2][1]= 0.084569; cbarray[28][3][1]= 1.451289;
cbarray[29][1][1]= 5.187717; cbarray[29][2][1]=-4.384724; cbarray[29][3][1]= 0.285782;
cbarray[30][1][1]= 6.601117; cbarray[30][2][1]= 0.734092; cbarray[30][3][1]=-1.451289;
cbarray[31][1][1]= 4.608140; cbarray[31][2][1]= 4.990316; cbarray[31][3][1]=-0.285782;
cbarray[32][1][1]= 2.408965; cbarray[32][2][1]= 6.308420; cbarray[32][3][1]=-0.787787;
cbarray[33][1][1]= 5.650418; cbarray[33][2][1]= 3.028036; cbarray[33][3][1]=-2.263549;
cbarray[34][1][1]= 0.947341; cbarray[34][2][1]= 5.846885; cbarray[34][3][1]=-3.337118;
cbarray[35][1][1]=-3.167641; cbarray[35][2][1]= 5.963667; cbarray[35][3][1]=-0.787787;
cbarray[36][1][1]=-5.187717; cbarray[36][2][1]= 4.384724; cbarray[36][3][1]=-0.285782;
cbarray[37][1][1]=-1.660299; cbarray[37][2][1]= 5.685678; cbarray[37][3][1]=-3.337118;
cbarray[38][1][1]=-6.641272; cbarray[38][2][1]=-0.084569; cbarray[38][3][1]=-1.451289;
cbarray[39][1][1]=-5.980361; cbarray[39][2][1]= 2.309007; cbarray[39][3][1]=-2.263549;
cbarray[40][1][1]=-4.643408; cbarray[40][2][1]=-4.419841; cbarray[40][3][1]=-2.263549;
cbarray[41][1][1]=-5.725840; cbarray[41][2][1]=-1.808038; cbarray[41][3][1]=-3.188360;
cbarray[42][1][1]=-1.085068; cbarray[42][2][1]=-6.552577; cbarray[42][3][1]=-1.451289;
cbarray[43][1][1]= 1.883898; cbarray[43][2][1]=-6.369032; cbarray[43][3][1]=-1.451289;
cbarray[44][1][1]=-2.466029; cbarray[44][2][1]=-5.385371; cbarray[44][3][1]=-3.337118;
cbarray[45][1][1]= 5.152449; cbarray[45][2][1]=-3.814249; cbarray[45][3][1]=-2.263549;
cbarray[46][1][1]= 5.904939; cbarray[46][2][1]=-1.089009; cbarray[46][3][1]=-3.188360;
cbarray[47][1][1]= 3.110577; cbarray[47][2][1]=-5.040619; cbarray[47][3][1]=-3.337118;
cbarray[48][1][1]= 4.366675; cbarray[48][2][1]= 2.622671; cbarray[48][3][1]=-4.502625;
cbarray[49][1][1]= 1.862773; cbarray[49][2][1]= 4.123417; cbarray[49][3][1]=-5.074189;
cbarray[50][1][1]= 4.523978; cbarray[50][2][1]= 0.078197; cbarray[50][3][1]=-5.074189;
cbarray[51][1][1]=-2.356477; cbarray[51][2][1]= 3.862577; cbarray[51][3][1]=-5.074189;
cbarray[52][1][1]=-4.656463; cbarray[52][2][1]= 2.064849; cbarray[52][3][1]=-4.502625;
cbarray[53][1][1]=-4.499161; cbarray[53][2][1]=-0.479624; cbarray[53][3][1]=-5.074189;
cbarray[54][1][1]=-0.179099; cbarray[54][2][1]= 2.897047; cbarray[54][3][1]=-6.147758;
cbarray[55][1][1]=-1.069376; cbarray[55][2][1]=-3.872887; cbarray[55][3][1]=-5.484256;
cbarray[56][1][1]= 1.538264; cbarray[56][2][1]=-3.711680; cbarray[56][3][1]=-5.484256;
cbarray[57][1][1]= 2.199175; cbarray[57][2][1]=-1.318104; cbarray[57][3][1]=-6.296516;
cbarray[58][1][1]=-2.020076; cbarray[58][2][1]=-1.578943; cbarray[58][3][1]=-6.296516;
cbarray[59][1][1]= 0.000000; cbarray[59][2][1]= 0.000000; cbarray[59][3][1]=-6.798521;
*/


cbarray[ 0][1][1]= 0.000000; cbarray[ 0][2][1]= 2.427051; cbarray[ 0][3][1]= 0.500000;
cbarray[ 1][1][1]=-0.809017; cbarray[ 1][2][1]= 2.118034; cbarray[ 1][3][1]= 1.000000;
cbarray[ 2][1][1]=-0.500000; cbarray[ 2][2][1]= 1.618034; cbarray[ 2][3][1]= 1.809017;
cbarray[ 3][1][1]= 0.809017; cbarray[ 3][2][1]= 2.118034; cbarray[ 3][3][1]= 1.000000;
cbarray[ 4][1][1]= 0.500000; cbarray[ 4][2][1]= 1.618034; cbarray[ 4][3][1]= 1.809017;
cbarray[ 5][1][1]= 1.618034; cbarray[ 5][2][1]= 1.809017; cbarray[ 5][3][1]= 0.500000;
cbarray[ 6][1][1]= 1.618034; cbarray[ 6][2][1]= 1.809017; cbarray[ 6][3][1]=-0.500000;
cbarray[ 7][1][1]= 0.000000; cbarray[ 7][2][1]= 2.427051; cbarray[ 7][3][1]=-0.500000;
cbarray[ 8][1][1]= 0.809017; cbarray[ 8][2][1]= 2.118034; cbarray[ 8][3][1]=-1.000000;
cbarray[ 9][1][1]=-0.809017; cbarray[ 9][2][1]= 2.118034; cbarray[ 9][3][1]=-1.000000;
cbarray[10][1][1]=-1.618034; cbarray[10][2][1]= 1.809017; cbarray[10][3][1]= 0.500000;
cbarray[11][1][1]=-1.618034; cbarray[11][2][1]= 1.809017; cbarray[11][3][1]=-0.500000;
cbarray[12][1][1]=-2.118034; cbarray[12][2][1]= 1.000000; cbarray[12][3][1]= 0.809017;
cbarray[13][1][1]=-1.000000; cbarray[13][2][1]= 0.809017; cbarray[13][3][1]= 2.118034;
cbarray[14][1][1]=-1.809017; cbarray[14][2][1]= 0.500000; cbarray[14][3][1]= 1.618034;
cbarray[15][1][1]=-0.500000; cbarray[15][2][1]= 0.000000; cbarray[15][3][1]= 2.427051;
cbarray[16][1][1]= 1.000000; cbarray[16][2][1]= 0.809017; cbarray[16][3][1]= 2.118034;
cbarray[17][1][1]= 0.500000; cbarray[17][2][1]= 0.000000; cbarray[17][3][1]= 2.427051;
cbarray[18][1][1]= 2.118034; cbarray[18][2][1]= 1.000000; cbarray[18][3][1]= 0.809017;
cbarray[19][1][1]= 1.809017; cbarray[19][2][1]= 0.500000; cbarray[19][3][1]= 1.618034;
cbarray[20][1][1]= 2.118034; cbarray[20][2][1]= 1.000000; cbarray[20][3][1]=-0.809017;
cbarray[21][1][1]= 2.427051; cbarray[21][2][1]= 0.500000; cbarray[21][3][1]= 0.000000;
cbarray[22][1][1]= 1.809017; cbarray[22][2][1]= 0.500000; cbarray[22][3][1]=-1.618034;
cbarray[23][1][1]= 0.500000; cbarray[23][2][1]= 1.618034; cbarray[23][3][1]=-1.809017;
cbarray[24][1][1]=-0.500000; cbarray[24][2][1]= 1.618034; cbarray[24][3][1]=-1.809017;
cbarray[25][1][1]= 1.000000; cbarray[25][2][1]= 0.809017; cbarray[25][3][1]=-2.118034;
cbarray[26][1][1]=-1.000000; cbarray[26][2][1]= 0.809017; cbarray[26][3][1]=-2.118034;
cbarray[27][1][1]=-2.118034; cbarray[27][2][1]= 1.000000; cbarray[27][3][1]=-0.809017;
cbarray[28][1][1]=-2.427051; cbarray[28][2][1]= 0.500000; cbarray[28][3][1]= 0.000000;
cbarray[29][1][1]=-1.809017; cbarray[29][2][1]= 0.500000; cbarray[29][3][1]=-1.618034;
cbarray[30][1][1]=-2.427051; cbarray[30][2][1]=-0.500000; cbarray[30][3][1]= 0.000000;
cbarray[31][1][1]=-1.809017; cbarray[31][2][1]=-0.500000; cbarray[31][3][1]= 1.618034;
cbarray[32][1][1]=-1.000000; cbarray[32][2][1]=-0.809017; cbarray[32][3][1]= 2.118034;
cbarray[33][1][1]=-2.118034; cbarray[33][2][1]=-1.000000; cbarray[33][3][1]= 0.809017;
cbarray[34][1][1]=-0.500000; cbarray[34][2][1]=-1.618034; cbarray[34][3][1]= 1.809017;
cbarray[35][1][1]= 1.000000; cbarray[35][2][1]=-0.809017; cbarray[35][3][1]= 2.118034;
cbarray[36][1][1]= 1.809017; cbarray[36][2][1]=-0.500000; cbarray[36][3][1]= 1.618034;
cbarray[37][1][1]= 0.500000; cbarray[37][2][1]=-1.618034; cbarray[37][3][1]= 1.809017;
cbarray[38][1][1]= 2.427051; cbarray[38][2][1]=-0.500000; cbarray[38][3][1]= 0.000000;
cbarray[39][1][1]= 2.118034; cbarray[39][2][1]=-1.000000; cbarray[39][3][1]= 0.809017;
cbarray[40][1][1]= 1.809017; cbarray[40][2][1]=-0.500000; cbarray[40][3][1]=-1.618034;
cbarray[41][1][1]= 2.118034; cbarray[41][2][1]=-1.000000; cbarray[41][3][1]=-0.809017;
cbarray[42][1][1]= 0.500000; cbarray[42][2][1]= 0.000000; cbarray[42][3][1]=-2.427051;
cbarray[43][1][1]=-0.500000; cbarray[43][2][1]= 0.000000; cbarray[43][3][1]=-2.427051;
cbarray[44][1][1]= 1.000000; cbarray[44][2][1]=-0.809017; cbarray[44][3][1]=-2.118034;
cbarray[45][1][1]=-1.809017; cbarray[45][2][1]=-0.500000; cbarray[45][3][1]=-1.618034;
cbarray[46][1][1]=-2.118034; cbarray[46][2][1]=-1.000000; cbarray[46][3][1]=-0.809017;
cbarray[47][1][1]=-1.000000; cbarray[47][2][1]=-0.809017; cbarray[47][3][1]=-2.118034;
cbarray[48][1][1]=-1.618034; cbarray[48][2][1]=-1.809017; cbarray[48][3][1]= 0.500000;
cbarray[49][1][1]=-0.809017; cbarray[49][2][1]=-2.118034; cbarray[49][3][1]= 1.000000;
cbarray[50][1][1]=-1.618034; cbarray[50][2][1]=-1.809017; cbarray[50][3][1]=-0.500000;
cbarray[51][1][1]= 0.809017; cbarray[51][2][1]=-2.118034; cbarray[51][3][1]= 1.000000;
cbarray[52][1][1]= 1.618034; cbarray[52][2][1]=-1.809017; cbarray[52][3][1]= 0.500000;
cbarray[53][1][1]= 1.618034; cbarray[53][2][1]=-1.809017; cbarray[53][3][1]=-0.500000;
cbarray[54][1][1]= 0.000000; cbarray[54][2][1]=-2.427051; cbarray[54][3][1]= 0.500000;
cbarray[55][1][1]= 0.500000; cbarray[55][2][1]=-1.618034; cbarray[55][3][1]=-1.809017;
cbarray[56][1][1]=-0.500000; cbarray[56][2][1]=-1.618034; cbarray[56][3][1]=-1.809017;
cbarray[57][1][1]=-0.809017; cbarray[57][2][1]=-2.118034; cbarray[57][3][1]=-1.000000;
cbarray[58][1][1]= 0.809017; cbarray[58][2][1]=-2.118034; cbarray[58][3][1]=-1.000000;
cbarray[59][1][1]= 0.000000; cbarray[59][2][1]=-2.427051; cbarray[59][3][1]=-0.500000;

          ADJ(0 , 1) = ADJ(1 , 0) = ADJ(0 , 3) = ADJ(3 , 0) =
          ADJ(0 , 7) = ADJ(7 , 0) = ADJ(1 , 2) = ADJ(2 , 1) =
          ADJ(1 ,10) = ADJ(10, 1) = ADJ(2 , 4) = ADJ(4 , 2) =
          ADJ(2 ,13) = ADJ(13, 2) = ADJ(3 , 4) = ADJ(4 , 3) =
          ADJ(3 , 5) = ADJ(5 , 3) = ADJ(4 ,16) = ADJ(16, 4) =
          ADJ(5 , 6) = ADJ(6  ,5) = ADJ(5 ,18) = ADJ(18, 5) =
          ADJ(6 , 8) = ADJ(8 , 6) = ADJ(6 ,20) = ADJ(20, 6) =
          ADJ(7 , 8) = ADJ(8 , 7) = ADJ(7 , 9) = ADJ(9 , 7) =
          ADJ(8 ,23) = ADJ(23, 8) = ADJ(9 ,11) = ADJ(11, 9) =
          ADJ(9 ,24) = ADJ(24, 9) = ADJ(10,11) = ADJ(11,10) =
          ADJ(10,12) = ADJ(12,10) = ADJ(11,27) = ADJ(27,11) =
          ADJ(12,14) = ADJ(14,12) = ADJ(12,28) = ADJ(28,12) =
          ADJ(13,14) = ADJ(14,13) = ADJ(13,15) = ADJ(15,13) =
          ADJ(14,31) = ADJ(31,14) = ADJ(15,17) = ADJ(17,15) =
          ADJ(15,32) = ADJ(32,15) = ADJ(16,17) = ADJ(17,16) =
          ADJ(16,19) = ADJ(19,16) = ADJ(17,35) = ADJ(35,17) =
          ADJ(18,19) = ADJ(19,18) = ADJ(18,21) = ADJ(21,18) =
          ADJ(19,36) = ADJ(36,19) = ADJ(20,21) = ADJ(21,20) =
          ADJ(20,22) = ADJ(22,20) = ADJ(21,38) = ADJ(38,21) =
          ADJ(22,25) = ADJ(25,22) = ADJ(22,40) = ADJ(40,22) =
          ADJ(23,24) = ADJ(24,23) = ADJ(23,25) = ADJ(25,23) =
          ADJ(24,26) = ADJ(26,24) = ADJ(25,42) = ADJ(42,25) =
          ADJ(26,29) = ADJ(29,26) = ADJ(26,43) = ADJ(43,26) =
          ADJ(27,28) = ADJ(28,27) = ADJ(27,29) = ADJ(29,27) =
          ADJ(28,30) = ADJ(30,28) = ADJ(29,45) = ADJ(45,29) =
          ADJ(30,33) = ADJ(33,30) = ADJ(30,46) = ADJ(46,30) =
          ADJ(31,32) = ADJ(32,31) = ADJ(31,33) = ADJ(33,31) =
          ADJ(32,34) = ADJ(34,32) = ADJ(33,48) = ADJ(48,33) =
          ADJ(34,37) = ADJ(37,34) = ADJ(34,49) = ADJ(49,34) =
          ADJ(35,36) = ADJ(36,35) = ADJ(35,37) = ADJ(37,35) =
          ADJ(36,39) = ADJ(39,36) = ADJ(37,51) = ADJ(51,37) =
          ADJ(38,39) = ADJ(39,38) = ADJ(38,41) = ADJ(41,38) =
          ADJ(39,52) = ADJ(52,39) = ADJ(40,41) = ADJ(41,40) =
          ADJ(40,44) = ADJ(44,40) = ADJ(41,53) = ADJ(53,41) =
          ADJ(42,43) = ADJ(43,42) = ADJ(42,44) = ADJ(44,42) =
          ADJ(43,47) = ADJ(47,43) = ADJ(44,55) = ADJ(55,44) =
          ADJ(45,46) = ADJ(46,45) = ADJ(45,47) = ADJ(47,45) =
          ADJ(46,50) = ADJ(50,46) = ADJ(47,56) = ADJ(56,47) =
          ADJ(48,49) = ADJ(49,48) = ADJ(48,50) = ADJ(50,48) =
          ADJ(49,54) = ADJ(54,49) = ADJ(50,57) = ADJ(57,50) =
          ADJ(51,52) = ADJ(52,51) = ADJ(51,54) = ADJ(54,51) =
          ADJ(52,53) = ADJ(53,52) = ADJ(53,58) = ADJ(58,53) =
          ADJ(54,59) = ADJ(59,54) = ADJ(55,56) = ADJ(56,55) =
          ADJ(55,58) = ADJ(58,55) = ADJ(56,57) = ADJ(57,56) =
          ADJ(57,59) = ADJ(59,57) = ADJ(58,59) = ADJ(59,58) = 1;

       break;


       case AB60:  /*   THE BALL WITH CENTRE   */
                   /*   --------------------   */


cbarray[ 0][1][1]= 0.000000; cbarray[ 0][2][1]= 0.000000; cbarray[ 0][3][1]= 6.798521;
cbarray[ 1][1][1]= 2.020076; cbarray[ 1][2][1]= 1.578943; cbarray[ 1][3][1]= 6.296516;
cbarray[ 2][1][1]= 1.069376; cbarray[ 2][2][1]= 3.872887; cbarray[ 2][3][1]= 5.484256;
cbarray[ 3][1][1]=-2.199175; cbarray[ 3][2][1]= 1.318104; cbarray[ 3][3][1]= 6.296516;
cbarray[ 4][1][1]=-1.538264; cbarray[ 4][2][1]= 3.711680; cbarray[ 4][3][1]= 5.484256;
cbarray[ 5][1][1]=-4.523978; cbarray[ 5][2][1]=-0.078197; cbarray[ 5][3][1]= 5.074189;
cbarray[ 6][1][1]=-4.366675; cbarray[ 6][2][1]=-2.622671; cbarray[ 6][3][1]= 4.502625;
cbarray[ 7][1][1]= 0.179099; cbarray[ 7][2][1]=-2.897047; cbarray[ 7][3][1]= 6.147758;
cbarray[ 8][1][1]=-1.862773; cbarray[ 8][2][1]=-4.123417; cbarray[ 8][3][1]= 5.074189;
cbarray[ 9][1][1]= 2.356477; cbarray[ 9][2][1]=-3.862577; cbarray[ 9][3][1]= 5.074189;
cbarray[10][1][1]= 4.499161; cbarray[10][2][1]= 0.479624; cbarray[10][3][1]= 5.074189;
cbarray[11][1][1]= 4.656463; cbarray[11][2][1]=-2.064849; cbarray[11][3][1]= 4.502625;
cbarray[12][1][1]= 5.725840; cbarray[12][2][1]= 1.808038; cbarray[12][3][1]= 3.188360;
cbarray[13][1][1]= 2.466029; cbarray[13][2][1]= 5.385371; cbarray[13][3][1]= 3.337118;
cbarray[14][1][1]= 4.643408; cbarray[14][2][1]= 4.419841; cbarray[14][3][1]= 2.263549;
cbarray[15][1][1]= 1.085068; cbarray[15][2][1]= 6.552577; cbarray[15][3][1]= 1.451289;
cbarray[16][1][1]=-3.110577; cbarray[16][2][1]= 5.040619; cbarray[16][3][1]= 3.337118;
cbarray[17][1][1]=-1.883898; cbarray[17][2][1]= 6.369032; cbarray[17][3][1]= 1.451289;
cbarray[18][1][1]=-5.904939; cbarray[18][2][1]= 1.089009; cbarray[18][3][1]= 3.188360;
cbarray[19][1][1]=-5.152449; cbarray[19][2][1]= 3.814249; cbarray[19][3][1]= 2.263549;
cbarray[20][1][1]=-5.650418; cbarray[20][2][1]=-3.028036; cbarray[20][3][1]= 2.263549;
cbarray[21][1][1]=-6.601117; cbarray[21][2][1]=-0.734092; cbarray[21][3][1]= 1.451289;
cbarray[22][1][1]=-4.608140; cbarray[22][2][1]=-4.990316; cbarray[22][3][1]= 0.285782;
cbarray[23][1][1]=-0.947341; cbarray[23][2][1]=-5.846885; cbarray[23][3][1]= 3.337118;
cbarray[24][1][1]= 1.660299; cbarray[24][2][1]=-5.685678; cbarray[24][3][1]= 3.337118;
cbarray[25][1][1]=-2.408965; cbarray[25][2][1]=-6.308420; cbarray[25][3][1]= 0.787787;
cbarray[26][1][1]= 3.167641; cbarray[26][2][1]=-5.963667; cbarray[26][3][1]= 0.787787;
cbarray[27][1][1]= 5.980361; cbarray[27][2][1]=-2.309007; cbarray[27][3][1]= 2.263549;
cbarray[28][1][1]= 6.641272; cbarray[28][2][1]= 0.084569; cbarray[28][3][1]= 1.451289;
cbarray[29][1][1]= 5.187717; cbarray[29][2][1]=-4.384724; cbarray[29][3][1]= 0.285782;
cbarray[30][1][1]= 6.601117; cbarray[30][2][1]= 0.734092; cbarray[30][3][1]=-1.451289;
cbarray[31][1][1]= 4.608140; cbarray[31][2][1]= 4.990316; cbarray[31][3][1]=-0.285782;
cbarray[32][1][1]= 2.408965; cbarray[32][2][1]= 6.308420; cbarray[32][3][1]=-0.787787;
cbarray[33][1][1]= 5.650418; cbarray[33][2][1]= 3.028036; cbarray[33][3][1]=-2.263549;
cbarray[34][1][1]= 0.947341; cbarray[34][2][1]= 5.846885; cbarray[34][3][1]=-3.337118;
cbarray[35][1][1]=-3.167641; cbarray[35][2][1]= 5.963667; cbarray[35][3][1]=-0.787787;
cbarray[36][1][1]=-5.187717; cbarray[36][2][1]= 4.384724; cbarray[36][3][1]=-0.285782;
cbarray[37][1][1]=-1.660299; cbarray[37][2][1]= 5.685678; cbarray[37][3][1]=-3.337118;
cbarray[38][1][1]=-6.641272; cbarray[38][2][1]=-0.084569; cbarray[38][3][1]=-1.451289;
cbarray[39][1][1]=-5.980361; cbarray[39][2][1]= 2.309007; cbarray[39][3][1]=-2.263549;
cbarray[40][1][1]=-4.643408; cbarray[40][2][1]=-4.419841; cbarray[40][3][1]=-2.263549;
cbarray[41][1][1]=-5.725840; cbarray[41][2][1]=-1.808038; cbarray[41][3][1]=-3.188360;
cbarray[42][1][1]=-1.085068; cbarray[42][2][1]=-6.552577; cbarray[42][3][1]=-1.451289;
cbarray[43][1][1]= 1.883898; cbarray[43][2][1]=-6.369032; cbarray[43][3][1]=-1.451289;
cbarray[44][1][1]=-2.466029; cbarray[44][2][1]=-5.385371; cbarray[44][3][1]=-3.337118;
cbarray[45][1][1]= 5.152449; cbarray[45][2][1]=-3.814249; cbarray[45][3][1]=-2.263549;
cbarray[46][1][1]= 5.904939; cbarray[46][2][1]=-1.089009; cbarray[46][3][1]=-3.188360;
cbarray[47][1][1]= 3.110577; cbarray[47][2][1]=-5.040619; cbarray[47][3][1]=-3.337118;
cbarray[48][1][1]= 4.366675; cbarray[48][2][1]= 2.622671; cbarray[48][3][1]=-4.502625;
cbarray[49][1][1]= 1.862773; cbarray[49][2][1]= 4.123417; cbarray[49][3][1]=-5.074189;
cbarray[50][1][1]= 4.523978; cbarray[50][2][1]= 0.078197; cbarray[50][3][1]=-5.074189;
cbarray[51][1][1]=-2.356477; cbarray[51][2][1]= 3.862577; cbarray[51][3][1]=-5.074189;
cbarray[52][1][1]=-4.656463; cbarray[52][2][1]= 2.064849; cbarray[52][3][1]=-4.502625;
cbarray[53][1][1]=-4.499161; cbarray[53][2][1]=-0.479624; cbarray[53][3][1]=-5.074189;
cbarray[54][1][1]=-0.179099; cbarray[54][2][1]= 2.897047; cbarray[54][3][1]=-6.147758;
cbarray[55][1][1]=-1.069376; cbarray[55][2][1]=-3.872887; cbarray[55][3][1]=-5.484256;
cbarray[56][1][1]= 1.538264; cbarray[56][2][1]=-3.711680; cbarray[56][3][1]=-5.484256;
cbarray[57][1][1]= 2.199175; cbarray[57][2][1]=-1.318104; cbarray[57][3][1]=-6.296516;
cbarray[58][1][1]=-2.020076; cbarray[58][2][1]=-1.578943; cbarray[58][3][1]=-6.296516;
cbarray[59][1][1]= 0.000000; cbarray[59][2][1]= 0.000000; cbarray[59][3][1]=-6.798521;
cbarray[60][1][1]= 0.0     ; cbarray[60][2][1]= 0.0     ; cbarray[60][3][1]= 0.0     ;

          ADJ(0 , 1) = ADJ(1 , 0) = ADJ(0 , 3) = ADJ(3 , 0) =
          ADJ(0 , 7) = ADJ(7 , 0) = ADJ(1 , 2) = ADJ(2 , 1) =
          ADJ(1 ,10) = ADJ(10, 1) = ADJ(2 , 4) = ADJ(4 , 2) =
          ADJ(2 ,13) = ADJ(13, 2) = ADJ(3 , 4) = ADJ(4 , 3) =
          ADJ(3 , 5) = ADJ(5 , 3) = ADJ(4 ,16) = ADJ(16, 4) =
          ADJ(5 , 6) = ADJ(6  ,5) = ADJ(5 ,18) = ADJ(18, 5) =
          ADJ(6 , 8) = ADJ(8 , 6) = ADJ(6 ,20) = ADJ(20, 6) =
          ADJ(7 , 8) = ADJ(8 , 7) = ADJ(7 , 9) = ADJ(9 , 7) =
          ADJ(8 ,23) = ADJ(23, 8) = ADJ(9 ,11) = ADJ(11, 9) =
          ADJ(9 ,24) = ADJ(24, 9) = ADJ(10,11) = ADJ(11,10) =
          ADJ(10,12) = ADJ(12,10) = ADJ(11,27) = ADJ(27,11) =
          ADJ(12,14) = ADJ(14,12) = ADJ(12,28) = ADJ(28,12) =
          ADJ(13,14) = ADJ(14,13) = ADJ(13,15) = ADJ(15,13) =
          ADJ(14,31) = ADJ(31,14) = ADJ(15,17) = ADJ(17,15) =
          ADJ(15,32) = ADJ(32,15) = ADJ(16,17) = ADJ(17,16) =
          ADJ(16,19) = ADJ(19,16) = ADJ(17,35) = ADJ(35,17) =
          ADJ(18,19) = ADJ(19,18) = ADJ(18,21) = ADJ(21,18) =
          ADJ(19,36) = ADJ(36,19) = ADJ(20,21) = ADJ(21,20) =
          ADJ(20,22) = ADJ(22,20) = ADJ(21,38) = ADJ(38,21) =
          ADJ(22,25) = ADJ(25,22) = ADJ(22,40) = ADJ(40,22) =
          ADJ(23,24) = ADJ(24,23) = ADJ(23,25) = ADJ(25,23) =
          ADJ(24,26) = ADJ(26,24) = ADJ(25,42) = ADJ(42,25) =
          ADJ(26,29) = ADJ(29,26) = ADJ(26,43) = ADJ(43,26) =
          ADJ(27,28) = ADJ(28,27) = ADJ(27,29) = ADJ(29,27) =
          ADJ(28,30) = ADJ(30,28) = ADJ(29,45) = ADJ(45,29) =
          ADJ(30,33) = ADJ(33,30) = ADJ(30,46) = ADJ(46,30) =
          ADJ(31,32) = ADJ(32,31) = ADJ(31,33) = ADJ(33,31) =
          ADJ(32,34) = ADJ(34,32) = ADJ(33,48) = ADJ(48,33) =
          ADJ(34,37) = ADJ(37,34) = ADJ(34,49) = ADJ(49,34) =
          ADJ(35,36) = ADJ(36,35) = ADJ(35,37) = ADJ(37,35) =
          ADJ(36,39) = ADJ(39,36) = ADJ(37,51) = ADJ(51,37) =
          ADJ(38,39) = ADJ(39,38) = ADJ(38,41) = ADJ(41,38) =
          ADJ(39,52) = ADJ(52,39) = ADJ(40,41) = ADJ(41,40) =
          ADJ(40,44) = ADJ(44,40) = ADJ(41,53) = ADJ(53,41) =
          ADJ(42,43) = ADJ(43,42) = ADJ(42,44) = ADJ(44,42) =
          ADJ(43,47) = ADJ(47,43) = ADJ(44,55) = ADJ(55,44) =
          ADJ(45,46) = ADJ(46,45) = ADJ(45,47) = ADJ(47,45) =
          ADJ(46,50) = ADJ(50,46) = ADJ(47,56) = ADJ(56,47) =
          ADJ(48,49) = ADJ(49,48) = ADJ(48,50) = ADJ(50,48) =
          ADJ(49,54) = ADJ(54,49) = ADJ(50,57) = ADJ(57,50) =
          ADJ(51,52) = ADJ(52,51) = ADJ(51,54) = ADJ(54,51) =
          ADJ(52,53) = ADJ(53,52) = ADJ(53,58) = ADJ(58,53) =
          ADJ(54,59) = ADJ(59,54) = ADJ(55,56) = ADJ(56,55) =
          ADJ(55,58) = ADJ(58,55) = ADJ(56,57) = ADJ(57,56) =
          ADJ(57,59) = ADJ(59,57) = ADJ(58,59) = ADJ(59,58) = 1;

       break;


       case SQ5:  /*   THE SQAURE WITH CENTRE   */
                  /*   ----------------------   */

          a    = 1.0/sqrt(2.0);
          cbarray[0][1][1]= a; cbarray[0][2][1]= a; cbarray[0][3][1]=0.;
          cbarray[1][1][1]= a; cbarray[1][2][1]=-a; cbarray[1][3][1]=0.;
          cbarray[2][1][1]=-a; cbarray[2][2][1]=-a; cbarray[2][3][1]=0.;
          cbarray[3][1][1]=-a; cbarray[3][2][1]= a; cbarray[3][3][1]=-0.000001;
          cbarray[4][1][1]=0.; cbarray[4][2][1]=0.; cbarray[4][3][1]= 0.000001;

          ADJ(0, 1) = ADJ(1 ,0) = ADJ(0 , 3) = ADJ(3 , 0) =
          ADJ(0, 4) = ADJ(4 ,0) = ADJ(1 , 2) = ADJ(2 , 1) =
          ADJ(1, 4) = ADJ(4 ,1) = ADJ(2 , 3) = ADJ(3 , 2) =
          ADJ(2, 4) = ADJ(4 ,2) = ADJ(3 , 4) = ADJ(4 , 3) = 1;

       break;

       default:
       printf("ERR* Can't calculate symmetry measure of %i polyhedron *ERR\n",
       number_case);
       exit(1);
       break;
       }
}




/**********************************************************************/



double determinant(matrix R)
{
    double tmp;
    tmp = 0.0;
    tmp += R[1][1]*(R[2][2]*R[3][3]-R[2][3]*R[3][2]);
    tmp += R[1][2]*(R[2][3]*R[3][1]-R[2][1]*R[3][3]);
    tmp += R[1][3]*(R[2][1]*R[3][2]-R[2][2]*R[3][1]);
    return(tmp);
}



/**********************************************************************/


void test_symmetry()
{
    int i, j, k;
    long ii, nn, mod, rem;
    // a unused
    //double a=2.0/sqrt(3.0);
    matrix T1;
    permuter *p = createPermuter(count, count);   

    rotation = alloc_matrix(3,3);

    while (nextPermutation(p)) {
	for (i = 0; i < count; i++) {
		best_permute[i] = p->_index[i];
	} 	
	eval_cube_symm();	
    }
/*
    nn = 1;
    for(i=0; i<count; i++)
       nn=nn*count;
    nn -= 1;
	
    for(ii=0; ii<nn; ii++)
       {
       mod = ii;
       for(i=0; i<count; i++)
          {
          rem = mod%count;
          mod = mod/count;
          best_permute[i] = rem;	   
          }	
       k = 0;
       for(i=0; i<count; i++)
          for(j=0; j<count; j++)
             if((best_permute[i] == best_permute[j]) && (i != j))
                { k=1; break; }
       if(k == 0) { 
          eval_cube_symm();
	}
       }
*/

    for(i=0; i<count; i++)
       {
       T1 = mul_mtx(rotation,3,3,cbarray[best_perm[i]],3,1);
       parray[i][1][1] = T1[1][1];
       parray[i][2][1] = T1[2][1];
       parray[i][3][1] = T1[3][1];
       free_matrix(T1);
       }

//    for(i=0; i<count; i++)
//       for(j=0; j<count; j++)
//         ADJ_OLD(i,j) = ADJ(best_perm[i],best_perm[j]);


    for(i=0; i<count; i++)
       {
       parray[i][1][1] *= norma;
       parray[i][2][1] *= norma;
       parray[i][3][1] *= norma;
       }
}

void test_symmetry_perms()
{
    int i, j, k, perms;
    long ii, nn, mod, rem;
    // a unused
    //double a=2.0/sqrt(3.0);
    matrix T1;     

    rotation = alloc_matrix(3,3);

    for (perms = 0; perms < numPermutations; perms++){		
	for (i = 0; i < count; i++) {
		best_permute[i] = permutations[perms][i];				
	} 					
	eval_cube_symm();	
    }

    for(i=0; i<count; i++)
       {
       T1 = mul_mtx(rotation,3,3,cbarray[best_perm[i]],3,1);
       parray[i][1][1] = T1[1][1];
       parray[i][2][1] = T1[2][1];
       parray[i][3][1] = T1[3][1];
       free_matrix(T1);
       }

    for(i=0; i<count; i++)
       {
       parray[i][1][1] *= norma;
       parray[i][2][1] *= norma;
       parray[i][3][1] *= norma;
       }
}



/**********************************************************************/



void test_symmetry_new()
{
    int i, j;
    matrix T1;

    rotation = alloc_matrix(3,3);

    for(i=0; i<count; i++)
       best_permute[i] = i;
    eval_cube_symm();

    for(i=0; i<count; i++)
       {	
       T1 = mul_mtx(rotation,3,3,cbarray[best_perm[i]],3,1);
       parray[i][1][1] = T1[1][1];
       parray[i][2][1] = T1[2][1];
       parray[i][3][1] = T1[3][1];	
       free_matrix(T1);
       }

    for(i=0; i<count; i++)
       {
       parray[i][1][1] *= norma;
       parray[i][2][1] *= norma;
       parray[i][3][1] *= norma;
       }

//    for(i=0; i<count; i++)
//       for(j=0; j<count; j++)
//         ADJ_OLD(i,j) = ADJ(i,j);

    free_matrix(rotation);
}



/**********************************************************************/




void eval_cube_symm()
{
    int i,j;
    double ss, sum_t, sum_pt, coef;
    matrix q, w, v, ut, R;
    matrix T1, T2;	    

    q = alloc_matrix(3,3);
    for(i=1; i<=3; i++)
        for(j=1; j<=3; j++)
            q[i][j]=0.0;

    for(i=0; i<count; i++)
       {
       T1 = transpose_mtx(parray[i],3,1);		
       T2 = mul_mtx(cbarray[best_permute[i]],3,1,T1,1,3);
       in_place_add_mtx(q,3,3,T2,3,3);
       free_matrix(T1);
       free_matrix(T2);
       }    



    w = alloc_matrix(3,1);
    v = alloc_matrix(3,3);
    if (!svdcmp(q,3,3,w,v)) {	
	free_matrix(q);
    	free_matrix(w);    
    	free_matrix(v);    
	return;
    }

    
    ut = transpose_mtx(q,3,3);
    R = mul_mtx(v,3,3,ut,3,3);

    
    sum_t  = 0.0;
    sum_pt = 0.0;
    for(i=0; i<count; i++)
       {
       T1  = mul_mtx(R,3,3,cbarray[best_permute[i]],3,1);

       sum_t  +=  (SQR(T1[1][1])+
                   SQR(T1[2][1])+
                   SQR(T1[3][1]));

       sum_pt += (parray[i][1][1]*T1[1][1]+
                  parray[i][2][1]*T1[2][1]+
                  parray[i][3][1]*T1[3][1]);
	
       free_matrix(T1);
       }

    sum_t  /= (double)(count);
    sum_pt /= (double)(count);
    coef    = sum_pt/sum_t;

    ss   = 0.0;
    for(i=0; i<count; i++)
       {
       T1  = mul_mtx(R,3,3,cbarray[best_permute[i]],3,1);

       ss += (SQR(parray[i][1][1]-T1[1][1]*coef) +
              SQR(parray[i][2][1]-T1[2][1]*coef) +
              SQR(parray[i][3][1]-T1[3][1]*coef));

       free_matrix(T1);
       }
    ss /= (double)(count);    

#ifdef DEBUG
   printf("%4.2f\n",ss);
   printf("%4.2f\n",determinant(R));
#endif

    if((ss < s) && (allowImproperRotations || determinant(R) > 0.0))
       {		

       s     = ss;
       norma = coef;
       for(i=1; i<=3; i++) {
          for(j=1; j<=3; j++) {
             rotation[i][j] = R[i][j];		
	   }	   
	}	
       for(i=0; i<count; i++)
          best_perm[i] = best_permute[i];
       }
	
    free_matrix(q);
    free_matrix(w);
    free_matrix(ut);
    free_matrix(v);
    free_matrix(R);
}


/**********************************************************************/


void copy_points_to_draw()
{
   int i;
   for(i=0; i<count; i++)
      {
      pdraw[i][0] = parray[i][1][1];
      pdraw[i][1] = parray[i][2][1];
      pdraw[i][2] = parray[i][3][1];
      }
}


/******************************************************************/


void copy_points_to_new_draw()
{
   int i;

   for(i=0; i<count; i++)
      {
      pndraw[i][0] = parray[i][1][1];
      pndraw[i][1] = parray[i][2][1];
      pndraw[i][2] = parray[i][3][1];
      }
}


/******************************************************************/


void denormalize_array(int n)
{
    int i;

    for(i=0;i<n;i++)
       {
       parray[i][1][1] = parray[i][1][1] * norm + x_avg;
       parray[i][2][1] = parray[i][2][1] * norm + y_avg;
       parray[i][3][1] = parray[i][3][1] * norm + z_avg;
      }
}


/**************************************************************************/


void print_output(int number_case, double sigma_input,
                  double h_b_ratio)
{
    int i;
    double ee;
    double x[MAXPOINTS],y[MAXPOINTS],z[MAXPOINTS];

    ee         = sigma_input/norm;
    shift_symm = 3.0*ee*ee;
    rms_symm   = ee*sqrt((12.0*ee*ee+8.0*s)/count);

    s *= 100;

    switch(number_case)
       {
       case SHAPE:
	printf("SHAPE = %7lf\n",s);
       break;
       case B4:
       printf("TETRAHEDRON = %7lf\n",s);
       break;

       case AB4:
       printf("TETRAHEDRON = %7lf\n",s);
       break;

       case B5:
       printf("BIPYRAMIDE = %7lf\n",s);
       break;

       case B5_:
       printf("EQUIDISTANT BIPYRAMIDE = %7lf\n",s);
       break;

       case AB5:
       printf("BIPYRAMIDE = %7lf\n",s);
       break;

       case AB5_:
       printf("EQUIDISTANT BIPYRAMIDE = %7lf\n",s);
       break;

       case PR:
       printf("TRIGONAL PRISM = %7lf\n",s);
       break;

       case APR:
       printf("TRIGONAL PRISM = %7lf\n",s);
       break;

       case D3H:
       printf("D3h SYMMETRY VALUE = %7lf\n", s);
       printf("D3h SYMMETRY h/b-RATIO %7lf\n", h_b_ratio);
       break;

       case AD3H:
       printf("D3h SYMMETRY VALUE = %7lf\n",s);
       printf("D3h SYMMETRY, h/b-RATIO %7lf\n", h_b_ratio);
       break;

       case PR_EQ:
       printf("TRIGONAL PRISM = %7lf\n",s);
       break;

       case APR_EQ:
       printf("TRIGONAL PRISM = %7lf\n",s);
       break;

       case B6:
       printf("OCTAHEDRON = %7lf\n",s);
       break;

       case AB6:
       printf("OCTAHEDRON = %7lf\n",s);
       break;

       case B8:
       printf("CUBE = %7lf\n",s);
       break;

       case AB8:
       printf("CUBE = %7lf\n",s);
       break;

       case B12:
       printf("ICOSAHEDRON = %7lf\n",s);
       break;

       case AB12:
       printf("ICOSAHEDRON = %7lf\n",s);
       break;

       case B20:
       printf("DODECAHEDRON = %7lf\n",s);
       break;

       case AB20:
       printf("DODECAHEDRON = %7lf\n",s);
       break;

       case B60:
       printf("THE BALL = %7lf\n",s);
       break;

       case AB60:
       printf("THE BALL = %7lf\n",s);
       break;

       case SQ5:
       printf("THE SQUARE = %7lf\n",s);
       break;
       }

    //printf("\nSYMMETRY BIAS:       %10lf\n",shift_symm);
    //printf("SYMMETRY RMS:        %10lf\n",rms_symm);
}


/******************************************************************/
/*   Minimization using Newton method  */
/***************************************/

#define MAXCYCLE  50
#define ERROR     0.00001
#define MINSYM    1.0e-14

void d3h_symmetry(int permut_flag, int n, int number_case, double* h_b_ratio)
{
    int i, j;
    double delta=0.0, dz=0.0, s_0=0.0, s_plus=0.0, s_minus=0.0, s_0_old=0.0;
    matrix pparray[MAXPOINTS];

    for (i=0; i< n; i++)
       pparray[i] = alloc_matrix(3,1);

    normalize_array(n);

    for(i=0;i<n;i++)
       {
       pparray[i][1][1] = parray[i][1][1];
       pparray[i][2][1] = parray[i][2][1];
       pparray[i][3][1] = parray[i][3][1];
       }

    delta   = 0.0001;
    dz      = 0.0;
    s_0_old = 1.0;

/**********************/

    for(j=0; j<MAXCYCLE; j++)
       {
       create_hedron_array(n, number_case);
       for(i=0; i<n/2; i++)
          cbarray[i][3][1] += dz;
       for(i=n/2; i<n; i++)
          cbarray[i][3][1] -= dz;
       *h_b_ratio = sqrt(2.0)*cbarray[0][3][1];
       normalize_cbarray(n);
       for(i=0; i<n; i++)
          {
          parray[i][1][1] = pparray[i][1][1];
          parray[i][2][1] = pparray[i][2][1];
          parray[i][3][1] = pparray[i][3][1];
          }
       s_0_old = s_0;
       copy_points_to_draw();
       if(permut_flag == 0) {
          test_symmetry();
       } else if (permut_flag == 1) {
          test_symmetry_new(); 
	} else if (permut_flag == 2) {	
	  test_symmetry_perms();
	}
       s_0 = s;
       if(s_0 < MINSYM) break;
       if(fabs((s_0-s_0_old)/s_0_old) < ERROR) break;
       s   = MAXDOUBLE;

/**********************/

       create_hedron_array(n, number_case);
       for(i=0; i<n/2; i++)
          cbarray[i][3][1] += (dz+delta);
       for(i=n/2; i<6; i++)
          cbarray[i][3][1] -= (dz+delta);
       normalize_cbarray(n);
       for(i=0; i<n; i++)
          {
          parray[i][1][1] = pparray[i][1][1];
          parray[i][2][1] = pparray[i][2][1];
          parray[i][3][1] = pparray[i][3][1];
          }
       if(permut_flag == 0) {
          test_symmetry();
       } else if (permut_flag == 1) { 
          test_symmetry_new();
	} else if (permut_flag == 2) {	
	  test_symmetry_perms();
	}
       s_plus = s;
       s      = MAXDOUBLE;

/**********************/

       create_hedron_array(n, number_case);
       for(i=0; i<n/2; i++)
          cbarray[i][3][1] += (dz-delta);
       for(i=n/2; i<n; i++)
          cbarray[i][3][1] -= (dz-delta);
       normalize_cbarray(n);
       for(i=0; i<n; i++)
          {
          parray[i][1][1] = pparray[i][1][1];
          parray[i][2][1] = pparray[i][2][1];
          parray[i][3][1] = pparray[i][3][1];
          }
       if(permut_flag == 0) {
          test_symmetry();
       } else if (permut_flag == 1) { 
          test_symmetry_new();
	} else if (permut_flag == 2) {	
	  test_symmetry_perms();
	}
       s_minus = s;
       s       = MAXDOUBLE;

/**********************/

       dz     -= 0.5*delta*(s_plus-s_minus)/(s_plus+s_minus-2.0*s_0);
       }

    if(j == MAXCYCLE)
       printf("ERR* CAN'T REACH THE GIVEN ACCURACY *ERR");
       s = s_0;
}


/******************************************************************/


void arun3_chirality(int permut_flag, int number_case, double sigma_input)
{
    double h_b_ratio=0.0;
    switch(number_case)
       {
       case D3H:
	      d3h_symmetry(permut_flag, count, number_case, &h_b_ratio);
       break;

       case AD3H:
          d3h_symmetry(permut_flag, count, number_case, &h_b_ratio);
       break;

       default:	
	create_hedron_array(count, number_case);
          normalize_array(count);          
          normalize_cbarray(count);
          copy_points_to_draw();	
       if(permut_flag == 0) {
          test_symmetry();
       } else if (permut_flag == 1) { 
          test_symmetry_new();
	} else if (permut_flag == 2) {	
	  test_symmetry_perms();
	}
		
       break;
       }

    copy_points_to_new_draw();

/*    denormalize_array(count);*/

    print_output(number_case, sigma_input, h_b_ratio);
}




/*******************************************************************/


#define PROGNAME argv[0]


int main(int argc, char *argv[])
{
    double sigma_input = 0.019;
    int    permute_flag;
    int    number_case;
    char *file_in;
    char *shape_in;

    if ((argc < 3) || (argc>7))
    {
      printf("usage: %s input_file case_flag [ref_shape if case_flag is 0] permute_flag [perms_file if permute_flag is 2]  [allow_improper]\n",PROGNAME);
      printf("Input Parameters:\n");
      printf("input_file - The shape to compare to a reference\n");
      printf("case_flag - The shape to compare to (0 for a user defined reference shape)\n");
      printf("optional: ref_shape - The reference shaps - only if case_flag = 0\n");
      printf("permute_flag: \n\t0 for enumeration on all permutations\n\t1 for identity perm only");
      printf("\n\t2 for enumartion on permutations in a given permutation file\n");
      printf("optional: perms_file - A permutations file, only valid if permute_flag = 2.\n");
      printf("\tThe file is of the format: num_permutations + a list of permutations\n");
      printf("\tThe permutations are performed on the reference shape (ref_shape)\n");
      printf("optional: allow_improper - \n\t0 - doesn't allow improper rotations (default) \n");
      printf("\n\t1 - allow improper rotations\n");
      exit(-1);
    }	

    file_in = argv[1];		    
    sscanf(argv[2], "%i", &number_case);	
    if (number_case == SHAPE) {
	shape_in = argv[3];
    	sscanf(argv[4], "%i", &permute_flag);
	if (permute_flag == 2) { 
 		permsFile = argv[5];	
		if(argc == 7) sscanf(argv[6], "%d", &allowImproperRotations);
	} else {
	    	if(argc == 6) sscanf(argv[5], "%d", &allowImproperRotations);
	}
	
    } else {		
    	sscanf(argv[3], "%i", &permute_flag);
	if (permute_flag == 2) { 
		permsFile = argv[4];		
		if(argc == 6) sscanf(argv[5], "%d", &allowImproperRotations);
	} else {
	    	if(argc == 5) sscanf(argv[4], "%d", &allowImproperRotations);
	}
    }	

    if (!read_input(number_case, file_in, shape_in)) {
      fprintf(stderr, "ERR* failed reading data from %s *ERR\n",file_in); exit(-1); }
	
    arun3_chirality(permute_flag, number_case, sigma_input);
    return 0;
}
