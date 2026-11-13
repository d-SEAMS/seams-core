#ifndef BABEL_ADAPTER_H
#define BABEL_ADAPTER_H

// Include Open Babel classes for OBMol and OBConversion
#include <openbabel/mol.h>
#include "Molecule.h"

using namespace OpenBabel;

/** 
 * Create a molecule from an OBMol
 * 
 * @param obmol The OpenBabel Molecule
 * @param replaceSym Whether to ignore atom names or not
 */
 Molecule* babel2Mol(OBMol &obmol, int reaplaceSym);

/**
 * Updates the coordinates of the OpenBabel Molecule according to the Molecule data
 * 
 * @param obmol The OpenBable molecule
 * @param outAtoms The output atoms' coordinates
 */
 void updateCoordinates(OBMol& obmol, double **outAtoms);

/**
 * Read a molecule from a file file
 * 
 * @param filename The file name to read from
 * @param format   The file format in Open Babel terms
 * @return The OBMol read from the file
 */
OBMol readMolecule (char *filename, const char *format);

/**
 * Write a mol to file
 *  
 * @param mol The molecule
 * @param format The output format
 * @param file The file to write output to
 * @param filename The file's name (for extension-finding purpose)
 */
void writeMolecule(OBMol& mol, const char *format, FILE* file, char *filename);

#endif
