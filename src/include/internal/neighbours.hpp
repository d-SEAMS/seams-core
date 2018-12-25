#ifndef __NEIGHBOURS_H_
#define __NEIGHBOURS_H_

#include <mol_sys.hpp>
#include <generic.hpp>

namespace nneigh{

	struct Jatom
	{

		int index;  // Index
		double r; 		// Distance from iatom

	};

	struct JList
	{
		std::vector<Jatom> n;
		int nearest_neighbours=0;
	};

	struct NeighbourList
	{

		std::vector<JList> iVector; 	// Collection of points

	};

	// Inefficient O(n^2) implementation of neighbour lists
	molSys::PointCloud<molSys::Point<double>, double> neighList(double rcutoff, molSys::PointCloud<molSys::Point<double>, double> *yCloud, 
		int typeI, int typeJ);

	// Inefficient O(n^2) implementation of neighbour lists 
	// You can only use this for neighbour lists with one type
	molSys::PointCloud<molSys::Point<double>, double> neighList(double rcutoff, molSys::PointCloud<molSys::Point<double>, double> *yCloud, 
		int typeI);

	// Inefficient O(n^2) implementation of neighbour lists 
	// You can only use this for neighbour lists with one type
	molSys::PointCloud<molSys::Point<double>, double> halfNeighList(double rcutoff, molSys::PointCloud<molSys::Point<double>, double> *yCloud, 
		int typeI=1);

	// Comparator for std::sort 
	inline bool compareByLength(const Jatom &a, const Jatom &b)
	{
    	return a.r < b.r;
	}

	// Clear neighbour list for the i^th atom if it is already full
	molSys::PointCloud<molSys::Point<double>, double> clearNeighList(molSys::PointCloud<molSys::Point<double>, double> *yCloud, int iatom);

}

#endif // __NEIGHBOURS_H_
