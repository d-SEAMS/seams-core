//-----------------------------------------------------------------------------------
// d-SEAMS - Deferred Structural Elucidation Analysis for Molecular Simulations
// SPDX-License-Identifier: MIT
//-----------------------------------------------------------------------------------

// Compiled with -fvisibility=hidden, the same way libseams is. The exported
// functions return this DSO's view of the exception typeinfo.

#include <euclidean_mic.hpp>
#include <linkcell.hpp>
#include <pore.hpp>

#include <typeinfo>

#if defined(__GNUC__) || defined(__clang__)
#define SEAMS_PROBE_API __attribute__((visibility("default")))
#else
#define SEAMS_PROBE_API
#endif

extern "C" SEAMS_PROBE_API const void *seams_hidden_range_type() {
  return &typeid(gen::mic::RangeError);
}

extern "C" SEAMS_PROBE_API const void *seams_hidden_linkcell_type() {
  return &typeid(linkcell::Error);
}

extern "C" SEAMS_PROBE_API const void *seams_hidden_box_type() {
  return &typeid(pore::BoxError);
}
