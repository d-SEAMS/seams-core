#ifndef SEAMS_CAGE_CANON_H_
#define SEAMS_CAGE_CANON_H_

#include <string>
#include <vector>

// Nauty (McKay and Piperno, J. Symbolic Comput. 60, 94 (2014);
// 10.1016/j.jsc.2013.09.003) produces a canonical adjacency certificate
// for the undirected graph of a ring set. findHC/findDDC stay the
// enumerators; this is the label-independent signature.

namespace cage {

bool nautyAvailable();

// Hex encoding of the canonical adjacency matrix. Empty if nauty is
// off or the ring set has no vertices.
std::string canonicalCertificate(const std::vector<std::vector<int>> &rings);

// True when the ring graph is isomorphic to the hexagonal prism (HC).
bool isHexagonalPrism(const std::vector<std::vector<int>> &rings);

bool sameCertificate(const std::vector<std::vector<int>> &a,
                     const std::vector<std::vector<int>> &b);

} // namespace cage

#endif
