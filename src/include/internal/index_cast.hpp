//-----------------------------------------------------------------------------------
// d-SEAMS - Deferred Structural Elucidation Analysis for Molecular Simulations
// SPDX-License-Identifier: MIT
//-----------------------------------------------------------------------------------
#pragma once

#include <concepts>
#include <cstddef>
#include <stdexcept>
#include <utility>

namespace seams {

// Integer conversion for a value that must fit the target type: an index
// or count that cannot be negative, or a size that must fit an int. A
// value out of range throws instead of wrapping to a huge index or a
// negative count. Not for hot loops; those keep their own index type.
template <std::integral To, std::integral From> constexpr To checkedCast(From value) {
  if (!std::in_range<To>(value)) {
    throw std::out_of_range("integer value out of range for its target type");
  }
  return static_cast<To>(value);
}

template <std::integral From> constexpr std::size_t toIndex(From value) {
  return checkedCast<std::size_t>(value);
}

} // namespace seams
