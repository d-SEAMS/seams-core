//-----------------------------------------------------------------------------------
// d-SEAMS - Deferred Structural Elucidation Analysis for Molecular Simulations
// SPDX-License-Identifier: MIT
//-----------------------------------------------------------------------------------

// Key functions for the exception types defined in this repository.
// Each destructor is the first non-inline virtual, so this translation
// unit emits the only typeinfo and vtable. libyodaLib exports them:
// default visibility on ELF, dllexport when this file is built with
// SEAMS_EH_EXPORT. linkcell::Error's key function is not in this file.
// It lives in linkcell's shared library linkcell_eh, which libyodaLib
// links. libseams is compiled with hidden visibility and does not emit
// a second copy of RangeError or BoxError.

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
