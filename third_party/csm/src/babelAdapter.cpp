#include <iostream>
#include <sstream>

// Include Open Babel classes for OBMol and OBConversion
#include <openbabel/mol.h>
#include <openbabel/atom.h>
#include <openbabel/bond.h>
#include <openbabel/obconversion.h>
#include <stdio.h>
#include "babelAdapter.h"

using namespace OpenBabel;
using namespace std;
extern "C" {
#include "Molecule.h"
	Molecule * allocateMolecule(int size);
	void initSimilarity(Molecule *m,int depth);		
	void replaceSymbols(Molecule* m);
}

const char *(ELEMENT_TABLE[]) = {
	"H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne", 
	"Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar", "K", "Ca",
	"Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", 
	"Ga", "Ge", "As", "Se",	"Br","Kr", "Rb", "Sr", "Y", "Zr",
	"Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn",
	"Sb", "Te", "I", "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd",
	"Pn", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb",
	"Lu", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg",
	"Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra",	"Ac", "Th",
	"Pa", "U", "Np", "Pu", "An", "Cn", "Bk", "Cf", "Es","Fm",
	"Md", "No", "Lr", "Rf", "Ha"
};

#define N_ELEMS 105

/** 
 * Create a molecule from an OBMol
 * 
 * @param obmol The OpenBabel Molecule
 * @param replaceSym Whether to ignore atom names or not
 */
Molecule* babel2Mol(OBMol &obmol, int replaceSym) {
	int numAtoms = obmol.NumAtoms();
	int i,j;
	Molecule *mol = allocateMolecule(numAtoms);

	for (i = 0; i < numAtoms; i++) {		
		OBAtom* atom = obmol.GetAtom(i + 1);	
		if (atom->GetAtomicNum() <= N_ELEMS && atom->GetAtomicNum() > 0) {
			mol->_symbol[i] = strdup(ELEMENT_TABLE[atom->GetAtomicNum() - 1]);		
		} else {
			mol->_symbol[i] = strdup(atom->GetType());		
		}
		mol->_valency[i] = atom->GetExplicitValence();
		mol->_adjacent[i] = (int*)malloc(mol->_valency[i] * sizeof(int));
		j = 0;
		for (OBBondIterator itr = atom->BeginBonds(); itr != atom->EndBonds(); itr++) {
			// Check if it's the first or second atom that is the neighbour			
			if ((*itr)->GetBeginAtomIdx() == atom->GetIdx()) {
				mol->_adjacent[i][j] = (*itr)->GetEndAtomIdx() - 1;	
			} else {
				mol->_adjacent[i][j] = (*itr)->GetBeginAtomIdx() - 1;	
			}
			j++;
		}
		mol->_pos[i][0] = atom->GetX();
		mol->_pos[i][1] = atom->GetY();
		mol->_pos[i][2] = atom->GetZ();
	}

	if (replaceSym) replaceSymbols(mol);

//	initSimilarity(mol,DEPTH_ITERATIONS);

	return mol;
}

/**
 * Updates the coordinates of the OpenBabel Molecule according to the Molecule data
 * 
 * @param obmol The OpenBable molecule
 * @param outAtoms The output atoms' coordinates
 */
void updateCoordinates(OBMol& obmol, double **outAtoms) {
	int numAtoms = obmol.NumAtoms();
	int i;
	for (i = 0; i < numAtoms; i++) {
		OBAtom* atom = obmol.GetAtom(i + 1);
		atom->SetVector(outAtoms[i][0], outAtoms[i][1], outAtoms[i][2]);
	}
}

/**
 * Read a molecule from a file file
 * 
 * @param filename The file name to read from
 * @param format   The file format in Open Babel terms (if NULL is given, attempt to get format 
 * 	  from the suffix)
 * @return The OBMol read from the file
 */
OBMol readMolecule (char *filename, const char *format) {
	OBMol mol;
	OBConversion conv;		
	if (format == NULL) {
		OBFormat* f = OBConversion::FormatFromExt(filename);
		if (f == NULL) {	
			printf("Error discovering format from filename %s\n", filename);
			exit(1);
		}	
		format = OBConversion::FormatFromExt(filename)->GetType().name();	
		if (!conv.SetInFormat(OBConversion::FormatFromExt(filename))) {
			printf ("Error setting input format to %s\n", format);
			exit(1);
		}
	} else {
		if (!conv.SetInFormat(format)) {
			printf ("Error setting input format to %s\n", format);
			exit(1);
		}
	}
	if (!conv.ReadFile(&mol, filename)) {
		printf ("Error reading file %s using OpenBabel with format %s\n", filename, format);
		exit(1);		
	}		
	return mol;
}

/**
 * Write a mol to file
 *  
 * @param mol The molecule
 * @param format The output format
 * @param file The file to write output to
 * @param filename The file's name (for extension-finding purpose)
 */
void writeMolecule(OBMol& mol, const char *format, FILE* file, char *filename) {
	ostringstream os;
	OBConversion conv;
	if (format == NULL) {
		OBFormat* f = OBConversion::FormatFromExt(filename);
		if (f != NULL) {	
			printf("Error discovering format from filename %s\n", filename);
			exit(1);
		}
		format = OBConversion::FormatFromExt(filename)->GetType().name();
		if (!conv.SetOutFormat(OBConversion::FormatFromExt(filename))) {
			printf ("Error setting output format to %s\n", format);
			exit(1);
		}		
	} else {
		if (!conv.SetOutFormat(format)) {
			printf ("Error setting output format to %s\n", format);
			exit(1);
		}
	}
	if (!conv.Write(&mol, &os)) {
		printf ("Error writing data file using OpenBabel with format %s\n", format);
		exit(1);
	}
	fprintf(file,os.str().c_str());
}
