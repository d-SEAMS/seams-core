#ifndef SEAMS_DOMAIN_HPP_
#define SEAMS_DOMAIN_HPP_

#include <cstdint>
#include <vector>

#include <mol_sys.hpp>

#ifdef SEAMS_HAS_MPI
#include <mpi.h>
#endif

/** @file domain.hpp
 *  @brief Spatial domain decomposition of one periodic frame across ranks.
 *
 *  Cells of the frame are ordered along a Hilbert curve and cut into one
 *  contiguous run per rank of near equal atom count, as Gadget splits its
 *  Peano-Hilbert keys (Springel, MNRAS 364, 1105, 2005). A rank's halo is
 *  every cell within a given distance of the cells it owns. A kernel whose
 *  answer for an atom depends only on atoms within that distance can run on
 *  the owned and halo atoms alone and agree with the whole frame on the
 *  owned ones.
 */

namespace seams::domain {

using Cloud = molSys::PointCloud<molSys::Point<double>, double>;

//! One rank's share of a frame
struct Share {
  std::vector<int> local;  //!< owned and halo atoms by frame index, ascending
  std::vector<char> owned; //!< nonzero where local[k] belongs to this rank
};

//! Position of cell (x, y, z) along the Hilbert curve over a 2^bits grid
std::uint64_t hilbertKey(std::uint32_t x, std::uint32_t y, std::uint32_t z,
                         int bits);

//! The atoms rank @a rank of @a nRanks owns, and every atom within @a halo
//! of one of them. Without a periodic box rank 0 owns every atom.
Share decompose(const Cloud &cloud, int rank, int nRanks, double halo);

//! The listed atoms of @a cloud, in that order, in the same box
Cloud subCloud(const Cloud &cloud, const std::vector<int> &atoms);

//! Primitive rings of the cutoff graph up to @a maxDepth members whose
//! lowest-indexed member rank @a rank owns, by frame index, as
//! primitive::ringNetwork lists them for the graph with ascending rows. Each
//! rank needs only its share, and the ranks' rings partition the frame's.
std::vector<std::vector<int>> rings(const Cloud &cloud, double cutoff,
                                    int maxDepth, int rank, int nRanks);

#ifdef SEAMS_HAS_MPI
//! Every rank's rings on every rank of @a comm, in the order rings() lists
//! them for a single rank. Collective over @a comm. Each rank builds the
//! whole list, so this costs every rank the size of the result, which must
//! fit in an int count of integers; linkcell's threads come from
//! RAYON_NUM_THREADS, every core unless a rank sets it.
std::vector<std::vector<int>> gatherRings(const Cloud &cloud, double cutoff,
                                          int maxDepth, MPI_Comm comm);
#endif

} // namespace seams::domain

#endif // SEAMS_DOMAIN_HPP_
