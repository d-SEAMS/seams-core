/*
 * Author: shadi lahham
 *
 * general functions used for parsing input from a file stream
 * supports custom XYZ format and PDB format
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>  //for ispace
#include <string.h>
#include "parseFunctions.h"

/*
 * removes all whitespace and stops at first non-whitespace character
 */
void eatWhitespace(FILE *in){

    char c;

    c = getc(in);
    while ( isspace( (unsigned char)c ))
        c = getc(in);
    ungetc(c,in);

}

/*
 * Reads a string from file stream
 * Allocates memory for string and returns
 */
char* readString(FILE *in){

    int len;
    char c, buffer[STRING_BUFFER_SIZE], *string, *q, *r;

    eatWhitespace(in);

    // read till next whitespace
    len = 0;

    c = getc(in);
	while( (c != EOF) && (! ( isspace( (unsigned char)c ) ) ) ){
		buffer[len++] = c;
		c = getc(in);

		// stop on overflow
		if (len+1 >= STRING_BUFFER_SIZE){
            break;

        }
	}
    buffer[len++] = '\0';

	// allocate space
	string = (char *)malloc(len+1 * sizeof(char) ); // +1 ?

	//copy buffer
	q = buffer;
	r = string;
	while (( *r++ = *q++) != '\0');

	return string;
}


// ******************************************************************************
//
//    PDB parse functions
//
//    format as specified by PDB reference:
//    http://www.rcsb.org/pdb/file_formats/pdb/pdbguide2.2/guide2.2_frame.html
//
// ******************************************************************************

int isPDBAtom(char * line){
	return ( ( strncmp(line,"ATOM",4) == 0 ) || ( strncmp(line,"HETATM",6) == 0 ) );
}

int isPDBEndModel(char * line){
	return ( strncmp(line,"ENDMDL",6) == 0 );
}

int isPDBconnect(char * line){
	return ( strncmp(line,"CONECT",6) == 0 );
}

/*
 * Counts number of atoms in PDB file
 * if file has more than one model just count the first model's atoms
 * rests file to start
 */
int countAtomsPDB(FILE *in){

	char line[LINE_LENGTH];
	int atoms = 0;

	while (fgets(line,LINE_LENGTH,in) ){
		if ( isPDBAtom (line) )
			atoms++;
		if ( isPDBEndModel (line) )
			break;
	}

	rewind(in);
	return(atoms);
}

/*
 * reads an PDB Atom symbol and x,y,z positions
 */
int readAtomPDB(FILE *in,char** sym_ptr,double* pos){

	char line[LINE_LENGTH],buffer[8],*sym;

	// read until encounter an atom line
	while (1) {
		if (! fgets(line,LINE_LENGTH,in) )
			return FALSE;

		if (isPDBAtom(line))
			break;
	}

	// allocate space for and read symbol
	sym = (char *)malloc(3 * sizeof(char) );
	strncpy(sym,line+12,2);
	sym[2]='\0';
	*sym_ptr = sym;

	// read x,y,z
	strncpy(buffer,line+30,8);
	pos[0] = atof(buffer);

	strncpy(buffer,line+38,8);
	pos[1] = atof(buffer);

	strncpy(buffer,line+46,8);
	pos[2] = atof(buffer);

	return TRUE;
}

/*
 * reads one line of PDB CONECT and updates valency and neighbours
 */
int readConnectivityLine(char *line, int *valency, int *neighbours,int size,int *curAtom){

	int neighbour;

	// skip CONECT
	line += 6;

	// read current atom number
	if ( ( sscanf(line,"%d",curAtom) != 1) || (*curAtom > size) )
		return FALSE;
	(*curAtom)--;

	// skip curAtom
	line += 5;

	// read the following neighbours , if any
	while (1){

		// eatWhitespace
		while ( isspace( (unsigned char)line[0] ) && (line[0] != '\n' ) )
			line++;

		if (line[0] == '\n')
			break;

		// scan number and skip digits in line
		sscanf(line,"%d",&neighbour);
		while (isdigit( (unsigned char)line[0] ) )
			line++;

		if (neighbour > size)
			return FALSE;

		neighbours[*valency] = neighbour -1; // -1 adjustment; PDB starts from 1 not 0
		(*valency)++;
	}

	return TRUE;
}

/*
 * reads PDB connectivity information and returns valency and neighbour list
 */
int readConnectivityPDB(FILE *in,int *valency,int *neighbours,int size,int *curAtom){

	char line[LINE_LENGTH],linehead[LINE_HEAD_LENGTH];
	fpos_t filePos;

	*valency = 0;

	// read until encounter a connect line
	while (1) {
		if (! fgets(line,LINE_LENGTH,in) )
			return FALSE;

		if (isPDBconnect(line))
			break;
	}

	// parse connect line
	if (! readConnectivityLine(line,valency,neighbours,size,curAtom) )
		return FALSE;

	// save head of line
	strncpy(linehead,line,LINE_HEAD_LENGTH);

	// read subsquent lines and parse if line head is the same
	while(1){

		// mark pos
		if (! ( fgetpos(in,&filePos) == 0 ) )
			break;

		// get next line
		if (! fgets(line,LINE_LENGTH,in) )
			break;

		// if not same head restore old position and stop
		if ( strncmp(linehead,line,LINE_HEAD_LENGTH) != 0){

			if (! ( fsetpos(in,&filePos) == 0 ) )
				return FALSE;          // error resetting - fatal
			break;
		}

		// parse the line
		if (! readConnectivityLine(line,valency,neighbours,size,curAtom) )
			return FALSE;

	}

	return TRUE;

}
