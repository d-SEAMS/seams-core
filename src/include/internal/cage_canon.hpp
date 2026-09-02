#ifndef SEAMS_CAGE_CANON_H_
#define SEAMS_CAGE_CANON_H_

#include <string>
#include <utility>
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

// Certificate of an arbitrary graph on n local vertices given by its edge
// list, with vertex `root` in its own colour cell so that rooted
// neighbourhoods with different centres get different certificates
// (root < 0 colours nothing). Empty if nauty is off, n is zero, or n
// exceeds the dense-graph limit.
std::string canonicalCertificateRooted(int n,
                                       const std::vector<std::pair<int, int>> &edges,
                                       int root);

// The same with a colour per vertex (an integer class such as the atom
// type): vertices of different colours never map onto each other, and the
// root still sits in a cell of its own. `colours` may be empty.
std::string canonicalCertificateColoured(int n,
                                         const std::vector<std::pair<int, int>> &edges,
                                         const std::vector<int> &colours, int root);

} // namespace cage

#endif
