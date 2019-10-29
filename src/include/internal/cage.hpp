#ifndef __CAGE_H_
#define __CAGE_H_
#include <vector>

// Namespace for cages
namespace cage {

// Type of a cage (a group of rings)
enum cageType { HexC, DoubleDiaC, Mixed };

// Each DDC has one equatorial ring and 6 peripheral rings
// Each HC has two basal planes and 3 prismatic planes
// Struct for point type
struct Cage {
  cageType type;          // type of the cage : can be DDC or HC
  std::vector<int> rings; // coordinates
};

} // namespace cage

#endif // __CAGE_H_