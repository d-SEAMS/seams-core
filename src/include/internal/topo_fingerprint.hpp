//-----------------------------------------------------------------------------------
// d-SEAMS - Deferred Structural Elucidation Analysis for Molecular Simulations
// SPDX-License-Identifier: MIT
//-----------------------------------------------------------------------------------
#ifndef SEAMS_TOPO_FINGERPRINT_H_
#define SEAMS_TOPO_FINGERPRINT_H_

#include <cstdint>
#include <map>
#include <string>
#include <utility>
#include <vector>

/** @file topo_fingerprint.hpp
 *  @brief Label-independent keys for bonded topologies.
 *
 *  A local key names the isomorphism class of the rooted bonded graph
 *  within a number of hops of one atom: two atoms share a key exactly when
 *  their neighbourhoods are the same graph with the same centre, whatever
 *  the atom indices. A frame key names the multiset of local keys together
 *  with the primitive ring census, so two configurations share a frame key
 *  when every atom has a counterpart with the same environment. This is
 *  the topological classification k-ART uses to index its event catalogue
 *  (Trochet, Beland, Joly, Brommer and Mousseau, Phys. Rev. B 91, 224106
 *  (2015)), here on the graphs the ring and cage code already builds.
 *
 *  With nauty linked the local key is the canonical certificate of the
 *  rooted graph (the centre in its own colour cell), an exact invariant.
 *  Without it the key is a Weisfeiler-Lehman colour refinement hash, which
 *  separates every pair of graphs the refinement can distinguish and is
 *  what the frame key uses in either build so that keys compare across
 *  hosts. The method is recorded on every result.
 */
namespace topo {

/// Neighbour rows by index, each row leading with the atom itself, as
/// nneigh::neighbourListByIndex returns them.
using Rows = std::vector<std::vector<int>>;

struct LocalKey {
  std::string key;     ///< canonical certificate or refinement hash
  std::string method;  ///< "nauty" or "wl"
  int vertices = 0;    ///< atoms in the neighbourhood, centre included
  int edges = 0;       ///< bonds among them
};

struct FrameFingerprint {
  std::string key;                    ///< hash of the sorted local keys and the ring census
  std::string method;                 ///< method of the local keys
  std::vector<std::string> atomKeys;  ///< one local key per atom
  std::map<std::string, int> classes; ///< local key -> number of atoms carrying it
  std::vector<int> ringCensus;        ///< ringCensus[s] = primitive rings of size s
  int hops = 0;
};

/// Atoms within `hops` bonds of `atom`, the atom itself first, then in
/// breadth-first order.
std::vector<int> hopNeighbourhood(const Rows &rows, int atom, int hops);

/// Weisfeiler-Lehman refinement hash of a graph given by local adjacency
/// lists; `root` (or -1) starts in its own colour; `colours` (empty or one
/// integer per vertex, such as the atom type) seeds the initial colours.
/// `rounds` refinements.
std::uint64_t wlHash(const std::vector<std::vector<int>> &adjacency, int root,
                     int rounds, const std::vector<int> &colours = {});

/// Key of the rooted neighbourhood of `atom`. `colours` is empty or one
/// class per row (atom type, species); coloured vertices never match
/// vertices of another colour.
LocalKey localKey(const Rows &rows, int atom, int hops,
                  const std::vector<int> &colours = {});

/// Keys of every atom, their histogram, the ring census up to
/// `maxRingSize`, and the frame key. `colours` as in localKey.
FrameFingerprint fingerprint(const Rows &rows, int hops = 2, int maxRingSize = 7,
                             const std::vector<int> &colours = {});

/// Hex string of a 64-bit hash.
std::string hex(std::uint64_t value);

/// A dictionary from local keys to labels: the keys of reference
/// structures (a polymorph library), so any molecule whose rooted
/// neighbourhood matches a reference gets that reference's name. Keys are
/// only comparable when they come from the same method, hop count and
/// colouring, which the header records.
struct KeyLibrary {
  std::string method;  ///< "nauty" or "wl"
  int hops = 0;
  std::map<std::string, std::string> labelOf;  ///< key -> label
};

/// Add every distinct key of `fp` under `label`; a key already present
/// under other labels carries all of them, sorted and joined by '|'.
void addToLibrary(KeyLibrary &lib, const FrameFingerprint &fp, const std::string &label);

/// Text form: a header line `# method M hops H`, then `key label` lines.
std::string writeLibrary(const KeyLibrary &lib);
KeyLibrary readLibrary(const std::string &text);

struct LibraryMatch {
  std::vector<std::string> labels;      ///< per atom; "" when no key matches
  std::map<std::string, int> counts;    ///< label -> atoms, "" for unmatched
  int matched = 0;
};

/// Look every atom key of `fp` up in `lib`. Throws std::invalid_argument
/// when the methods or hop counts differ.
LibraryMatch matchLibrary(const FrameFingerprint &fp, const KeyLibrary &lib);

} // namespace topo

#endif // SEAMS_TOPO_FINGERPRINT_H_
