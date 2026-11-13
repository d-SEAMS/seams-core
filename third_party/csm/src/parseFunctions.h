/*
 * Author: shadi lahham
 *
 * general functions used for parsing input from a file stream
 *
 */

#ifndef PARSEFUNCTIONS_H
#define PARSEFUNCTIONS_H

#define STRING_BUFFER_SIZE 200  /* maximal length of atom symbol */
#define LINE_LENGTH 1000        /* length of a line of input in a PDB file */
#define LINE_HEAD_LENGTH 11     /* length of CONECT and first digit in a PDB file = 6+5*/
#define TRUE 1
#define FALSE 0

void eatWhitespace(FILE *in);

char* readString(FILE *in);

int countAtomsPDB(FILE *in);

int readAtomPDB(FILE *in,char** sym_ptr,double* pos);

int readConnectivityPDB(FILE *in,int *valency,int *neighbours,int size,int *curAtom);

#endif
