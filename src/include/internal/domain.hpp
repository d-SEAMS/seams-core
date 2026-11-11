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
 *  The box is cut into slabs, pencils, or blocks, whichever keeps the most
 *  owned volume for a halo of the given width. A rank owns one block. Its
 *  halo is every atom within that width of the block, periodic boundaries
 *  included, and nothing further. A kernel whose answer for an atom depends
 *  only on atoms within that width can run on the owned and halo atoms alone
 *  and agree with the whole frame on the owned ones.
 */

namespace seams::domain {

using Cloud = molSys::PointCloud<molSys::Point<double>, double>;

//! One rank's share of a frame
struct Share {
  std::vector<int> local;  //!< owned and halo atoms by frame index, ascending
  std::vector<char> owned; //!< nonzero where local[k] belongs to this rank
  std::vector<int> owner;  //!< rank owning local[k]
};

//! What the ring analyses read off a frame's primitive rings
struct RingSummary {
  //! Rings of each size, indexed by member count, 0 through maxDepth
  std::vector<long long> sizeCount;
  //! Largest ring size each atom belongs to, 1 when it is in none, the rule
  //! ring::assignPolygonType applies. Filled on the root rank only.
  std::vector<int> atomType;
};

//! The summary of @a ringList over a frame of @a nAtoms atoms
RingSummary summarize(const std::vector<std::vector<int>> &ringList,
                      int nAtoms, int maxDepth);

//! Position of cell (x, y, z) along the Hilbert curve over a 2^bits grid
std::uint64_t hilbertKey(std::uint32_t x, std::uint32_t y, std::uint32_t z,
                         int bits);

//! Graph hops a primitive ring of @a maxDepth members can reach. A hop is
//! at most one cutoff, so the halo is this times the cutoff.
inline int ringHops(int maxDepth) {
  const int maxLvl = maxDepth > 0 ? maxDepth / 2 : 0;
  if (maxLvl <= 0) {
    return 0;
  }
  return maxLvl + (maxLvl - 1) / 2;
}

//! The atoms rank @a rank of @a nRanks owns, and every atom within @a halo
//! of its block. Without a periodic box rank 0 owns every atom.
Share decompose(const Cloud &cloud, int rank, int nRanks, double halo);

//! Ranks that own a halo atom of @a share, sorted, without this rank.
//! The relation is symmetric: a cell within @a halo of both sides is
//! listed by each owner.
std::vector<int> haloRanks(const Share &share);

//! The listed atoms of @a cloud, in that order, in the same box
Cloud subCloud(const Cloud &cloud, const std::vector<int> &atoms);

//! Primitive rings of the cutoff graph up to @a maxDepth members whose
//! lowest-indexed member rank @a rank owns, by frame index, as
//! primitive::ringNetwork lists them for the graph with ascending rows. Each
//! rank needs only its share, and the ranks' rings partition the frame's.
std::vector<std::vector<int>> rings(const Cloud &cloud, double cutoff,
                                    int maxDepth, int rank, int nRanks);

//! Rings of a cloud that is already one share, in share order. @a owned
//! selects the sources. Member indices are local.
std::vector<std::vector<int>> ownedRings(const Cloud &local, double cutoff,
                                         int maxDepth,
                                         const std::vector<char> &owned);

#ifdef SEAMS_HAS_MPI
//! Owned atoms plus the halo received from the neighbouring blocks.
struct Work {
  Cloud cloud;          //!< share.local order: index k is share.local[k]
  Share share;
  double commMs = 0;    //!< ghost exchange, after a barrier
  double haloFraction = 0; //!< halo atoms / local atoms, on this rank
};

//! Build @a rank's working cloud. Ghost coordinates arrive from the owners
//! over the distributed graph; the owned coordinates stay local.
Work gatherWork(const Cloud &cloud, double halo, MPI_Comm comm);

//! Every rank's rings on every rank of @a comm, in the order rings() lists
//! them for a single rank. Collective over @a comm. Each rank builds the
//! whole list, so this costs every rank the size of the result, which must
//! fit in an int count of integers; linkcell's threads come from
//! RAYON_NUM_THREADS, every core unless a rank sets it.
std::vector<std::vector<int>> gatherRings(const Cloud &cloud, double cutoff,
                                          int maxDepth, MPI_Comm comm);

//! summarize() of the whole frame's rings, with nothing whole-frame sized
//! leaving a rank. Each rank searches its own domain and halo for the rings
//! whose lowest-indexed member it owns. Ring counts are summed with
//! MPI_Allreduce, so every rank gets them. A halo atom's type travels to its
//! owner over MPI_Alltoallv, and only owned atoms in some ring reach @a root
//! through MPI_Gatherv, so atomType is filled on @a root alone. Collective
//! over @a comm.
RingSummary reduceRings(const Cloud &cloud, double cutoff, int maxDepth,
                        MPI_Comm comm, int root = 0);
#endif

} // namespace seams::domain

#endif // SEAMS_DOMAIN_HPP_
