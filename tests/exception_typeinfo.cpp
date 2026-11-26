//-----------------------------------------------------------------------------------
// d-SEAMS - Deferred Structural Elucidation Analysis for Molecular Simulations
// SPDX-License-Identifier: MIT
//-----------------------------------------------------------------------------------

// libc++ on Apple compares std::type_info by address. libstdc++ compares the
// type name, so a catch across two DSOs can succeed here while the addresses
// still differ. This test compares the addresses, which is the check that
// fails when a hidden DSO emits its own typeinfo.

#include <euclidean_mic.hpp>
#include <linkcell.hpp>
#include <pore.hpp>

#include <cstdio>
#include <typeinfo>

extern "C" const void *seams_hidden_range_type();
extern "C" const void *seams_hidden_linkcell_type();
extern "C" const void *seams_hidden_box_type();

static int same(const char *name, const void *here, const void *hidden) {
  if (here == hidden) {
    return 0;
  }
  std::fprintf(stderr, "%s typeinfo differs across DSOs (%p vs %p)\n", name, here, hidden);
  return 1;
}

int main() {
  int failed = 0;
  failed |= same("gen::mic::RangeError", &typeid(gen::mic::RangeError), seams_hidden_range_type());
  failed |= same("linkcell::Error", &typeid(linkcell::Error), seams_hidden_linkcell_type());
  failed |= same("pore::BoxError", &typeid(pore::BoxError), seams_hidden_box_type());
  return failed == 0 ? 0 : 1;
}
