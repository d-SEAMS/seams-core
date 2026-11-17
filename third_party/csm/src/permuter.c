/*
 * Author: shadi lahham
 *
 * generates permutations in the range 1 .. size
 * the permutations generated are limited to those containing
 * closed permutation groups of size groupSize
 *
 * permuter can generate permutations from a given integer array
 *
 */

#include <stdlib.h>  // for NULL
#include <stdio.h>   // for printf
#include "permuter.h"


/*
 * creates a permuter of size 'size' (with groupSize permutation subgroups)
 */

permuter* createPermuter(int size,int groupSize){
    // allocate permuter
    permuter *p = (permuter *)malloc(sizeof(permuter));

    //set size
    p->_size = size;

    //set groups
    p->_groupSize = groupSize;

    //allocate index
    p->_index = (int*)malloc(size * sizeof(int));

#ifdef USE_ORDER2_PERMUTER
	if (groupSize != 2) {
		printf("Error - operation order must be 2");
		exit(1);
	}

	p->_stack = (StackData *)malloc((size + 1) * sizeof(StackData));
	p->_inStack = 0;
#else
	//allocate used
	p->_used = (int*)malloc(size * sizeof(int));
#endif

    //init index
    resetPermuter(p);
    p->_firstPermutation = TRUE;

    return p;
}

void resetPermuter(permuter *p){
	int i;
	for (i=0; i<p->_size; i++){
		p->_index[i] = i;
	}
	p->_firstPermutation = TRUE;	
}

/*
 * swaps a with b in int array swaps
 */
void swap(int* s, int a, int b){
	int temp=s[a];
	if (a == b) return;
	s[a]=s[b];
	s[b]=temp;
}

/*
 * generates the next permutation regardless of groupSize
 */
int permute(permuter *p){

	int* str = p->_index;
	int len = p->_size;
	int key=len-1;
	int newkey=len-1;

	/* The key value is the first value from the
	   end which is smaller than the value to
	   its immediate right */
	while( (key>0)&&(str[key] <= str[key-1])){
		key--;
	}
	key--;

	/*If key<0 the data is in reverse sorted order,
	which is the last permutation.  */
	if(key <0) return 0;

	/*str[key+1] is greater than str[key] because
	of how key was found. If no other is greater,
	str[key+1] is used*/

	newkey=len-1;
	while((newkey > key) && (str[newkey] <= str[key])){
		newkey--;
	}
	swap(str, key,newkey);

	/* variables len and key are used to walk through
	   the tail, exchanging pairs from both ends of
	   the tail.  len and key are reused to save memory*/
	len--;
	key++;

	/* The tail must end in sorted order to produce the
	  next permutation.*/
	while(len>key){
		swap(str,len,key);
		key++;
		len--;
	}
	return 1;
}

/*
 * get the next valid permutation
 * in case of error or reached end of permutations, returns 0 (FALSE)
 */
int nextPermutation(permuter *p){
	// first permuation is unique, no need to permute

	if (p->_firstPermutation) {
		p->_firstPermutation = FALSE;
	} else {
		if (! permute(p)) return FALSE;	
	}

	return TRUE;
}

/*
 * free memory of the permuter structure
 */
void freePermuter(permuter *p){

    	// free index
	free(p->_index);

	// free used
	free(p->_used);

	//free permuter
	free(p);
};
