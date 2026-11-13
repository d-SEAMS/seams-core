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

#ifndef PERMUTER_H
#define PERMUTER_H

#define TRUE 1
#define FALSE 0

typedef struct permuter {
    int _size;
    int _groupSize;
    int* _index;			
    int* _used;             // internally used for marking
    int _firstPermutation;  // boolean flag - the first permuation is unique, not permuted
} permuter;

permuter* createPermuter(int size,int groupSize);

int nextPermutation(permuter *p);

void resetPermuter(permuter *p);

void freePermuter(permuter *p);

#endif
