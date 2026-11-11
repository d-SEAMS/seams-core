/*
 * Author: shadi lahham
 *
 * A simple molecule structure with :
 *  - atom positions XYZ
 *  - atom symbols H|C|Na|Cl .. etc
 *  - adjacency matrix to represent connectivity 1|0
 *  - valency of each atom
 *  - similarity: d-similar atoms have the same symbol and same structure
 *                of children and grandchildren up to a depth of d
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h> //for strcmp,strlen
#include <ctype.h>  //for ispace
#include <math.h>   //for sqrt
#include "Molecule.h"
#include "parseFunctions.h"


// ************************************************************
//       function declarations
// ************************************************************
void initSimilarity(Molecule *m,int depth);
void setMarked(Molecule *m,int state);
int isSimilar(Molecule *m,int a,int b);


// ************************************************************
//       implementation
// ************************************************************

/*
 * allocates memory for the molecule structure
 */
Molecule * allocateMolecule(int size){
    int i;

    // allocate molecule
    Molecule *m = (Molecule *)malloc(sizeof(Molecule));

    //set size
    m->_size = size;

    //allocate positions
	m->_pos = (double **)malloc(size * sizeof(double*));
	for (i=0;i<size;i++)
		m->_pos[i] = (double *)malloc(DIM * sizeof(double));

    //allocate symbols
	m->_symbol = (char **)malloc(size * sizeof(char*));
	// individual strings allocated on reading

    //allocate adjacency
    m->_adjacent = (int **)malloc(size * sizeof(int*));
	// individual neighbours allocated on reading

    //allocate valency
    m->_valency = (int*)malloc(size * sizeof(int));
    for ( i=0;  i< size ;  i++ ){
		m->_valency[i] = 0;
	}

    //allocate similarity
    m->_similar = (int*)malloc(size * sizeof(int));

    //allocate marked array
    m->_marked = (int*)malloc(m->_size * sizeof(int));

    return m;
};

/*
 * replace atom symbols with 'XX' - for unknown
 */
void replaceSymbols(Molecule* m){

	int i;

	// set new symbols 'XX'
	char *sym;

    // free old symbols
    for (i=0;i<m->_size;i++){
        if (m->_symbol[i]){
    		free(m->_symbol[i]);
        }
	}

    for (i=0;i<m->_size;i++){
    	sym = (char *)malloc(3 * sizeof(char) );
    	sym[0]='X';
    	sym[1]='X';
		sym[2]='\0';
		m->_symbol[i] = sym;
	}

}

/*
 * creates a molecule from an input data file
 */
Molecule* createMolecule(FILE *in,FILE *err,int replaceSym){

	Molecule* m;
	int size,i;

	    // read connectivity
    int atomNum,neighbour;
    char c;

    // allocate temporary buffer
    int *buff;
    int pos;
	int valency = 0;
    int j;

    // read size
    if (fscanf(in,"%d",&size)!=1){
        fprintf(err,"Input Error: Number of atoms not supplied!\n");
        return NULL;
    }
    if (size == 0)
		return NULL;

    // allocate molecule
    m = allocateMolecule(size);

    // read atoms
    for (i=0; i<size; i++){

        m->_symbol[i] = readString(in); // allocates space for symbol

        if(fscanf(in,"%lf%lf%lf",&(m->_pos[i][0]),&(m->_pos[i][1]),&(m->_pos[i][2]))!=3){
            fprintf(err,"Input Error: Failed reading input for atom %d\n",i /* +0 */ +1);
            freeMolecule(m);
            return NULL;
        }

    }

    if (replaceSym)
    	replaceSymbols(m);

	buff = (int*)malloc(size * sizeof(int));

    for (i=0; i<size; i++){

        // read and check atom number
        fscanf(in,"%d",&atomNum);
        if (atomNum != i /* +0 */ +1){
            fprintf(err,"Input Error: Failed reading connectivity for atom %d\n",i /* +0 */ +1);
            freeMolecule(m);
            return NULL;
        }

        
        pos = 0;

        // read neighbour numbers till newline
        while(1){

            // eat whitespace except newline
            c = getc(in);
            while ( isspace( (unsigned char)c ) && (c != '\n' ) )
                c = getc(in);

            // termination
            if (c == '\n')
                break;

            // otherwise put back read char
            ungetc(c,in);

            // read neighbour number
            if( (fscanf(in,"%d",&neighbour) !=1) || (neighbour > /* >= */ size) ) {
                fprintf(err,"Input Error: Failed reading connectivity for atom %d\n",i /* +0 */ +1);
                freeMolecule(m);
                return NULL;
            }

            // write neighbour to buffer
            buff[pos++] = neighbour  /* -0 */ -1;

            // increment valancy
            valency++;

        }
        m->_valency[i] = valency;

        // allocate memory for neighbours and copy
        m->_adjacent[i] = (int*)malloc(valency * sizeof(int));

        for (j=0; j<valency; j++)
        	m->_adjacent[i][j] = buff[j];

    }

    free(buff);

    initSimilarity(m,DEPTH_ITERATIONS);

	

    return (m);
}

/*
 * creates a molecule from a PDB input data file
 */
Molecule* createMoleculePDB(FILE *in,FILE *err,int replaceSym){

    int size,i;
	int *neighbours;
    int valency,curAtom;

	Molecule* m;
	int j;	

    // read size
    size = countAtomsPDB(in);
    if (size == 0)
		return NULL;

    // allocate molecule
    m = allocateMolecule(size);

    // read atoms
    for (i=0; i<size; i++){

		// note: readAtom allocates space for symbol

		if (! readAtomPDB(in,&(m->_symbol[i]),m->_pos[i]) ){
            fprintf(err,"Input Error: Failed reading input for atom %d\n",i);
            freeMolecule(m);
            return NULL;
        }

    }

    if (replaceSym)
    	replaceSymbols(m);

    // read connectivity
    neighbours = (int*)malloc(size * sizeof(int));

    for (i=0; i<size; i++){

		if ( (! readConnectivityPDB(in,&valency,neighbours,size,&curAtom) ) ||
			(curAtom < 0) ){
            fprintf(err,"Input Error: Failed reading connectivity element number %d\n",i+1);
            freeMolecule(m);
            return NULL;
        }

        m->_valency[curAtom] = valency;

        // allocate memory for neighbours and copy
        m->_adjacent[curAtom] = (int*)malloc(valency * sizeof(int));
        
        for (j=0; j<valency; j++)
        	m->_adjacent[curAtom][j] = neighbours[j];

    }

    free(neighbours);

    initSimilarity(m,DEPTH_ITERATIONS);

    return (m);
}

/*
 * Creates a new Molecule from selected atoms of the source Molecule (src)
 *
 * ! assumes selectedAtoms are present in src
 * ! This implies that the selectedAtomsSize is smaller or equal to source size
 *
 */
Molecule* copyMolecule(Molecule *src, int* selectedAtoms, int selectedAtomsSize, int updateSimilarity ){

	int i,j,old_i,len;
	int *buff;
	int *inversSelected;
	int *newGroups;


	// allocate molecule
    Molecule* dest = allocateMolecule(selectedAtomsSize);
    if (!dest){
    	return NULL; // allocation failed
    }

	//// allocate temporary buffers
	buff = (int*)malloc(dest->_size * sizeof(int));
	inversSelected = (int*)malloc(src->_size * sizeof(int));
	newGroups = (int*)malloc(src->_groupNum * sizeof(int));

	// init inversSelected - is the inverse of selectedAtoms (index->element,element->index)
	for ( i=0;  i< src->_size ;  i++ ) {
		inversSelected[i] = -1;	
	}
	for (i = 0; i < src->_groupNum; i++) {
		newGroups[i] = 0;
	}
	for ( i=0;  i< dest->_size ;  i++ )
		inversSelected[selectedAtoms[i]] = i;

	dest->_groupNum = 0;

	// main loop
	for (i=0; i<dest->_size; i++){

		// update adjacency and valency
		int pos = -1;
		int valency = 0;


		old_i = selectedAtoms[i];

		// copy pos
		for (j=0; j<DIM; j++){
			dest->_pos[i][j] = src->_pos[old_i][j];
		}

		// copy similar
		if (newGroups[src->_similar[old_i] - 1] == 0) {
			dest->_groupNum++;
			newGroups[src->_similar[old_i] - 1] = dest->_groupNum;
		}
		dest->_similar[i] = newGroups[src->_similar[old_i] - 1];

		// copy symbol - allocate the same size as old
		len = strlen(src->_symbol[old_i]);
		dest->_symbol[i] = (char *)malloc((len+1) * sizeof(char) );
		strcpy(dest->_symbol[i],src->_symbol[old_i]);

    	// for each item in src adjacency list
		for ( j=0;  j< src->_valency[old_i] ;  j++ ){

			// if item in selectedAtoms add to buffer
			pos = inversSelected[src->_adjacent[old_i][j]];
			if (pos != -1)
				buff[valency++] = pos;
		}

		dest->_valency[i] = valency;

		// allocate memory for adjacent and copy buffer
        dest->_adjacent[i] = (int*)malloc(valency * sizeof(int));
        for (j=0; j<valency; j++)
        	dest->_adjacent[i][j] = buff[j];

	}

	free(inversSelected);
	free(buff);
	free(newGroups);

	if (updateSimilarity)
		initSimilarity(dest,DEPTH_ITERATIONS);

	return(dest);
}

/*
 * breaks down atoms into similar groups by symbol and graph structure
 * depth -  is the desired maximal depth to take into account
 */
void initSimilarity(Molecule *m,int depth){

    int i,j,k,groupNum;
	int *group;
    int *subGroup;


    groupNum = 1;
    setMarked(m,FALSE);

    // break into initial groups by symbol and valancy
    for (i=0; i<m->_size ; i++){

        if (m->_marked[i])
            continue;

        for (j=0; j<m->_size ; j++){
            if ( (m->_marked[j]) || (m->_valency[i] != m->_valency[j]) || (strcmp(m->_symbol[i],m->_symbol[j])!=0) )
                 continue;

             m->_similar[j] = groupNum;
             m->_marked[j] = TRUE;
        }
        groupNum++;
    }

    group = (int*)malloc(m->_size * sizeof(int));     //temporary buffer
    subGroup = (int*)malloc(m->_size * sizeof(int));  //temporary buffer

    // iteratively refine the breakdown into groups
    for (k=0; k < depth ; k++){
	    int updated = FALSE;
		int len = 0; //group size

    	setMarked(m,FALSE); // reset marked

    	for (i=0; i < m->_size; i++){
			int subLen = 0; //newGroup size

	    	if (m->_marked[i])
	    		continue;

	    	// mark self as done
	    	m->_marked[i] = TRUE;

	    	// create the items in the group of i
	    	for(j=0; j < m->_size; j++)
	    		if ((j!=i) && (m->_similar[j] == m->_similar[i])){
	    			group[len++] = j;
	    		}

	    	// for each item in the group check if it can be split or not
	    	for(j=0; j < len; j++){
	    		if (isSimilar(m,group[j],i))
	    			// mark similar item as done
	    			m->_marked[group[j]] = TRUE;
	    		else
	    			// add to new subGroup
	    			subGroup[subLen++] = group[j];
	    	}

	    	// give subGroup members a new id
	    	for(j=0; j < subLen; j++){
	    		updated = TRUE;
	    		m->_similar[subGroup[j]] = groupNum;
	    	}

	    	if (updated == TRUE)
	    		groupNum++;
		}
    }

    m->_groupNum = groupNum -1;

    free(group);
    free(subGroup);
}

/*
 * sets the general use marked array to state (TRUE|FALSE)
 */
void setMarked(Molecule *m,int state){
    int i;
    for (i=0; i<m->_size; i++)
        m->_marked[i] = state;
}

/*
 * Atom a is similar to Atom b if for each neighbour of i, j has a similar neighbour
 */
int isSimilar(Molecule *m,int a,int b){

    int i,j,found = TRUE;

	// alloc and init temporary buffer
    int *mark = (int*)malloc(m->_size * sizeof(int));
    for ( i=0;  i<m->_size;  i++ )
    	mark[i] = 0;

	// for each of i's neighbours
	for ( i=0;  i<m->_valency[a];  i++ ){

		found = FALSE;

		for ( j=0;  j<m->_valency[b];  j++ ){
			if (mark[j])
				continue;

			if (m->_similar[m->_adjacent[a][i]] == m->_similar[m->_adjacent[b][j]]){
				found = TRUE;
				mark[j] = TRUE;
				break;
			}

		}

		if (!found)
			break;

	}

	free(mark);
	return(found);
}

/*
 * returns the number of elements in group 'num'
 * if there is no such group returns zero
 */
int getGroupSize(Molecule *m,int num){
	int i,j;
	i = 0;
	for (j=0; j<m->_size ; j++){
		if (m->_similar[j] == num){
			i++;
		}
	}
	return i;
}

/*
 * retrieves the similarity group number 'num' into the supplied buffer
 * returns the number of elements in that group
 * if there is no such group returns zero (buffer irrelevant)
 */
int getGroup(Molecule *m,int num,int* buff){
	int i,j;
	i = 0;
	for (j=0; j<m->_size ; j++){
		if (m->_similar[j] == num){
			buff[i] = j;
			i++;
		}
	}
	return i;
}

/*
 * retrieves the size of the largest similarity group
 */
int getMaxGroupSize(Molecule *m){

	int i,max=0;

	// alloc and init counting buffer
    int *count = (int*)malloc(m->_groupNum * sizeof(int));
    for ( i=0;  i<m->_groupNum;  i++ )
    	count[i] = 0;

	// count the number of members of each group
	for ( i=0;  i< m->_groupNum ;  i++ ){
		count[m->_similar[i] - 1]++;
	}

	// get the maximum
	for ( i=0;  i< m->_groupNum ;  i++ ){
		max = count[i] > max ? count[i] : max;
	}

	free(count);
	return max;
}

/*
 * Creates a new Molecule from m by removing atoms who's symbol is in the
 * remove list
 */
Molecule* stripAtoms(Molecule *m, char** removeList, int removeListSize, int updateSimilarity){

	Molecule* newM;

	int i, j, count, hits;

    int *selected = (int*)malloc(m->_size * sizeof(int));

	// find atoms not in removeList
    count = 0;

	for ( i=0;  i< m->_size ;  i++ ){
		hits = 0;
		for ( j=0;  j< removeListSize ;  j++ ) {
			if( strcmp( m->_symbol[i], removeList[j] ) == 0 ) {
				hits++;
				break;
			}
		}
		if (hits == 0){
			selected[count] = i;
			count ++;
		}
	}

	// return a new Molecule copy
	newM = copyMolecule(m,selected,count,updateSimilarity);
	free(selected);
	return newM;
}

/*
 * Normalizes the position of atoms of the molecule
 * returns one [TRUE] if successful, zero[FALSE] otherwise
 */
int normalizeMolecule(Molecule *m){

	double tmp,x_avg, y_avg, z_avg,norm;
	int i;

	x_avg = y_avg = z_avg = 0.0;

	for(i=0; i< m->_size; i++){
		x_avg += m->_pos[i][0];
		y_avg += m->_pos[i][1];
		z_avg += m->_pos[i][2];
	}
	x_avg /= (double)(m->_size);
	y_avg /= (double)(m->_size);
	z_avg /= (double)(m->_size);

	norm = 0.0;
	for(i=0; i< m->_size; i++){
		tmp = SQR(m->_pos[i][0]-x_avg) +
		      SQR(m->_pos[i][1]-y_avg) +
		      SQR(m->_pos[i][2]-z_avg);
		norm += tmp;
	}
	// normalize to 1 and not molecule size
	//norm = sqrt(norm / (double)m->_size);
	norm = sqrt(norm);

	if(norm < MINDOOUBLE)
		return FALSE;

	for(i=0; i< m->_size; i++){
		m->_pos[i][0] = ((m->_pos[i][0] - x_avg) / norm);
		m->_pos[i][1] = ((m->_pos[i][1] - y_avg) / norm);
		m->_pos[i][2] = ((m->_pos[i][2] - z_avg) / norm);
	}

	return TRUE;
}

/*
 * prints the molecule
 */
void printMolecule(Molecule *m){
    int i,j;

    // print molecule
    printf("molecule:\n");

    // print size
    printf("size = %d\n",m->_size);

    // print symbols and positions
	for (i=0;i<m->_size;i++)
        printf("%s %lf %lf %lf\n",m->_symbol[i],m->_pos[i][0],m->_pos[i][1],m->_pos[i][2]);

    // print adjacency & valency
    printf("\nconnectivity:\n");
   	for (i=0;i<m->_size;i++){
        printf("%d ",i+1);
       	for (j=0;j<m->_valency[i];j++)
        	printf("%d ",m->_adjacent[i][j]+1);
        printf("\t\tTotal = %d\n",m->_valency[i]);
    }

    printf("\nFULL similar subtree:\n");
   	for (i=0;i<m->_size;i++){
        printf("%d ",i+1);
        for (j=0;j<m->_size;j++){
            if (m->_similar[i] == m->_similar[j])
                printf("%d ",j+1);
        }
        printf("\n");
    }

    printf("\nsimilar subtree:\n");
    setMarked(m,FALSE);
   	for (i=0;i<m->_size;i++){

        if (m->_marked[i])
            continue;

        for (j=0; j<m->_size ; j++){
            if ( (m->_marked[j]) || (m->_similar[i] != m->_similar[j]) )
                 continue;

             printf("%d ",j+1);
             m->_marked[j] = TRUE;
        }
        printf("\n");
    }

    printf("\nsimilar array:\n");
    for (i=0;i<m->_size;i++){
        printf("%2d|",i+1);
    }
    printf("\n");
    for (i=0;i<m->_size;i++){
        printf("%2d|",m->_similar[i]);
    }
    printf("\n");

};

/*
 * prints the molecule - short version
 */
void printMoleculeBasic(Molecule *m){
    int i,j;

    // print size
    printf("%d\n",m->_size);

    // print symbols and positions
	for (i=0;i<m->_size;i++)
        printf("%2s %10lf %10lf %10lf\n",m->_symbol[i],m->_pos[i][0],m->_pos[i][1],m->_pos[i][2]);

    // print adjacency
   	for (i=0;i<m->_size;i++){
        printf("%d ",i+1);
       	for (j=0;j<m->_valency[i];j++)
        	printf("%d ",m->_adjacent[i][j]+1);
        printf("\n");
    }

};

/*
 * prints only the similar section of the Molecule
 */
void printMoleculeSimilar(Molecule *m){
    int i,j;

    printf("similar subtree:\n");
    setMarked(m,FALSE);
   	for (i=0;i<m->_size;i++){

        if (m->_marked[i])
            continue;

        for (j=0; j<m->_size ; j++){
            if ( (m->_marked[j]) || (m->_similar[i] != m->_similar[j]) )
                 continue;

             printf("%d ",j+1);
             m->_marked[j] = TRUE;
        }
        printf("\n");
    }

    printf("\nsimilar array:\n");
    for (i=0;i<m->_size;i++){
        printf("%2d|",i+1);
    }
    printf("\n");
    for (i=0;i<m->_size;i++){
        printf("%2d|",m->_similar[i]);
    }
    printf("\n");

};

/*
 * prints the Molecule with detailed information for debugging
 */
void printMoleculeDebug(Molecule *m){
    int i,j;

	// print similarity
	printf("Equivalent Groups:\n");
    setMarked(m,FALSE);
   	for (i=0;i<m->_size;i++){

        if (m->_marked[i])
            continue;

        for (j=0; j<m->_size ; j++){
            if ( (m->_marked[j]) || (m->_similar[i] != m->_similar[j]) )
                 continue;

             printf("%d ",j+1);
             m->_marked[j] = TRUE;
        }
        printf("\n");
    }

	printf("\n========DEBUG INFORMATION========\n");

    // print molecule
    printf("molecule:\n");

    // print size
    printf("size = %d\n",m->_size);

    // print symbols and positions
	for (i=0;i<m->_size;i++)
        printf("%s %lf %lf %lf\n",m->_symbol[i],m->_pos[i][0],m->_pos[i][1],m->_pos[i][2]);

    // print adjacency & valency
    printf("\nconnectivity:\n");
   	for (i=0;i<m->_size;i++){
        printf("%d ",i+1);
       	for (j=0;j<m->_valency[i];j++)
        	printf("%d ",m->_adjacent[i][j]+1);
        printf("\t\tTotal = %d\n",m->_valency[i]);
    }

    printf("\nComplete Equivalent Groups for each atom:\n");
   	for (i=0;i<m->_size;i++){
        printf("%d ",i+1);
        for (j=0;j<m->_size;j++){
            if (m->_similar[i] == m->_similar[j])
                printf("%d ",j+1);
        }
        printf("\n");
    }

};

/*
 * for debug purposes ... prints only the similar section of the Molecule
 */
void printMoleculeDebug2(Molecule *m){
    int i,j;

    printf("Breakdown into groups:\n");
    setMarked(m,FALSE);
   	for (i=0;i<m->_size;i++){

        if (m->_marked[i])
            continue;

        for (j=0; j<m->_size ; j++){
            if ( (m->_marked[j]) || (m->_similar[i] != m->_similar[j]) )
                 continue;

             printf("%d ",j+1);
             m->_marked[j] = TRUE;
        }
        printf("\n");
    }
    printf("\n");

};

/*
 * free memory of the molecule structure
 */
void freeMolecule(Molecule *m){
    int i;

    // free positions
    for (i=0;i<m->_size;i++){
		free(m->_pos[i]);
	}
	free(m->_pos);

    // free symbols
    for (i=0;i<m->_size;i++){
        if (m->_symbol[i]){
    		free(m->_symbol[i]);
        }
	}
	free(m->_symbol);

    // free adjacency
    for (i=0;i<m->_size;i++){
    	if (m->_adjacent[i]){
			free(m->_adjacent[i]);
		}
	}
	free(m->_adjacent);

    // free valency
	free(m->_valency);

    // free similar
	free(m->_similar);

    // free marked
	free(m->_marked);

	//free molecule
	free(m);

};
