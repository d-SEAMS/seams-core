#ifndef _GEOMETRY_H
#define _GEOMETRY_H

#include "molecular_system.h"

class CLayer {
	public:
		CLayer();
		virtual ~CLayer();

		// Check whether the given coordinate (r) is within the layer (midpoint r_layer)
		bool atomInsideLayer(double r, double r_layer, double dr);

};

// Class for volume 
class CVolume {
	public:
		CVolume();
		virtual ~CVolume();

		// Volume limits (may or may not be set)
		double xlo, xhi, ylo, yhi, zlo, zhi;
		// This is true if the volume has not been defined by the user (default)
		bool notset;
		// Flag to check if the atom is inside the given volume or not 
		int volFlag;
		// Total number of atoms of Itype
        int n_iatoms;
        // Total number of atoms of Jtype
        int n_jatoms;

        // Dynamically allocated array for indices with particles of type I and J
        int* iIndex;
        int* jIndex;

		// Check whether the atom is inside the user-defined volume
		bool atomInsideVol(class CMolecularSystem& molSys, int iatom, double xlo = 0.0, double xhi = 0.0, double ylo=0.0, double yhi=0.0, double zlo=0.0, double zhi=0.0);
		bool atomInsideVol(class CMolecularSystem& molSys, int iatom);

		// Check the volume limits
		void checkVolume();
		// Check if the atom is within a particular dimension range
		void atomCoordLimits(double r_atom, double r_min, double r_max);

		// Get list of atoms of a particular type in a user-defined volume
		void getAtomListI(class CMolecularSystem& molSys, int typeI=-1, double xlo = 0.0, double xhi = 0.0, double ylo=0.0, double yhi=0.0, double zlo=0.0, double zhi=0.0);
		void getAtomListIJ(class CMolecularSystem& molSys, int typeI=-1, int typeJ=-1, double xlo = 0.0, double xhi = 0.0, double ylo=0.0, double yhi=0.0, double zlo=0.0, double zhi=0.0);

};

class CGeneric {
	public:
		CGeneric();
		virtual ~CGeneric();
    	// Get absolute relative distance from wrapped coordinates
		double getAbsDistance(int, int, class CMolecularSystem& molSys);
		double getAbsDistance(int iatom, class CMolecularSystem* frameOne, class CMolecularSystem* frameTwo);
};

#endif
