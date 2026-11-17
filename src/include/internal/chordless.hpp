#ifndef SEAMS_CHORDLESS_H_
#define SEAMS_CHORDLESS_H_

#include <cstddef>
#include <vector>

namespace chordless {

//! Largest cycle the meet-in-the-middle search takes. A half is at most
//! six steps.
constexpr int kMaxSize = 12;

/**
 * @brief Chordless cycles of one size on a CSR graph.
 *
 * A cycle is chordless when no edge joins two of its vertices that are
 * not neighbours on the cycle. Each cycle is listed once, starting at its
 * lowest vertex, in the direction whose second vertex is the smaller of
 * that vertex's two cycle neighbours. The list is sorted
 * lexicographically, so it does not depend on the thread count.
 *
 * @param offsets n + 1 row offsets.
 * @param cols Column indices. Self entries must not appear.
 * @param size Cycle length, from 3 through kMaxSize.
 * @param rings When not null, receives count * size vertex indices.
 * @param per_atom When not null, n counts of cycles through each vertex.
 * @return The number of cycles.
 */
std::size_t cycles(const std::vector<int> &offsets, const std::vector<int> &cols, int size,
                   std::vector<int> *rings, int *per_atom);

} // namespace chordless

#endif
