//-----------------------------------------------------------------------------------
// d-SEAMS - Deferred Structural Elucidation Analysis for Molecular Simulations
// SPDX-License-Identifier: MIT
//-----------------------------------------------------------------------------------

// Key functions for the exception types defined in this repository.
// Each destructor is the first non-inline virtual, so this translation
// unit emits the only typeinfo and vtable, and libyodaLib exports them.
// linkcell::Error's key function is defined in that library; libyodaLib
// exports the typeinfo because it links the object. libseams is compiled
// with hidden visibility and does not emit a second copy.

#include <euclidean_mic.hpp>
#include <pore.hpp>

namespace gen {
namespace mic {
RangeError::~RangeError() = default;
}
} // namespace gen

namespace pore {
BoxError::~BoxError() = default;
}
