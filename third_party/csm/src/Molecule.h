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

#ifndef MOLECULE_H
#define MOLECULE_H

#define DIM 3
#define MINDOOUBLE  1e-8
#define LINE_BUFFER_SIZE 1000  /* maximal length of line of input */
#define DEPTH_ITERATIONS 200   /* maximal depth to descend to when checking similarity */
#define TRUE 1
#define FALSE 0

#define SQR(x)      ((x) * (x))


typedef struct MoleculeTag {
    int _size;
    double** _pos;           // atom positions XYZ
    char** _symbol;          // atom symbols
    int** _adjacent;         // represent connectivity
    int* _valency;           // valency of each atom
    int* _similar;           // similarity
    int* _marked;            // for marking atoms - general use
    int  _groupNum;          // the number of groups of similarity
} Molecule;

Molecule* createMolecule(FILE *in,FILE *err,int replaceSym);

Molecule* createMoleculePDB(FILE *in,FILE *err,int replaceSym);

Molecule* copyMolecule(Molecule *src, int* selectedAtoms, int selectedAtomsSize, int updateSimilarity );

int getGroup(Molecule *m,int num,int* buff);

int getGroupSize(Molecule *m,int num);

int getMaxGroupSize(Molecule *m);

Molecule* stripAtoms(Molecule *m, char** removeList, int removeListSize, int updateSimilarity);

int normalizeMolecule(Molecule *m);

void printMolecule(Molecule *m);

void printMoleculeBasic(Molecule *m);

void printMoleculeSimilar(Molecule *m);

void printMoleculeDebug(Molecule *m);

void printMoleculeDebug2(Molecule *m);

void freeMolecule(Molecule *m);

#endif
