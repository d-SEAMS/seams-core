//-----------------------------------------------------------------------------------
// d-SEAMS - Deferred Structural Elucidation Analysis for Molecular Simulations
// SPDX-License-Identifier: MIT
//-----------------------------------------------------------------------------------

// Compiled with -fvisibility=hidden where the compiler has that flag, the
// same way libseams is. The exported functions are this DSO's view of the
// exception types. On Windows the probe is a DLL and must dllexport them;
// a visibility attribute is not an export.

#include <euclidean_mic.hpp>
#include <linkcell.hpp>
#include <pore.hpp>

#include <typeinfo>

#if defined(_WIN32)
#define SEAMS_PROBE_API __declspec(dllexport)
#elif defined(__GNUC__) || defined(__clang__)
#define SEAMS_PROBE_API __attribute__((visibility("default")))
#else
#define SEAMS_PROBE_API
#endif

extern "C" SEAMS_PROBE_API const void *seams_hidden_range_type() {
  return &typeid(gen::mic::RangeError);
}

extern "C" SEAMS_PROBE_API void seams_hidden_throw_range() {
  throw gen::mic::RangeError("cross-dso-range");
}

extern "C" SEAMS_PROBE_API const void *seams_hidden_linkcell_type() {
  return &typeid(linkcell::Error);
}

extern "C" SEAMS_PROBE_API void seams_hidden_throw_linkcell() {
  throw linkcell::Error("cross-dso-linkcell");
}

extern "C" SEAMS_PROBE_API const void *seams_hidden_box_type() {
  return &typeid(pore::BoxError);
}

extern "C" SEAMS_PROBE_API void seams_hidden_throw_box() {
  throw pore::BoxError("cross-dso-box");
}
