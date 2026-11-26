//-----------------------------------------------------------------------------------
// d-SEAMS - Deferred Structural Elucidation Analysis for Molecular Simulations
// SPDX-License-Identifier: MIT
//-----------------------------------------------------------------------------------

// libc++ on Apple compares std::type_info by address. libstdc++ compares the
// type name, so a catch across two DSOs can succeed here while the addresses
// still differ. The address check is the one that fails when a hidden DSO
// emits its own typeinfo. It is Itanium (Linux and macOS). MSVC gives each
// module its own type_info object and matches the decorated name, so the
// addresses differ when the catch still hits. That catch runs on every host,
// including Windows.

#include <euclidean_mic.hpp>
#include <linkcell.hpp>
#include <pore.hpp>

#include <cstdio>
#include <cstring>
#include <exception>
#include <typeinfo>

extern "C" const void *seams_hidden_range_type();
extern "C" void seams_hidden_throw_range();
extern "C" const void *seams_hidden_linkcell_type();
extern "C" void seams_hidden_throw_linkcell();
extern "C" const void *seams_hidden_box_type();
extern "C" void seams_hidden_throw_box();

static int same(const char *name, const void *here, const void *hidden) {
  if (here == hidden) {
    return 0;
  }
  std::fprintf(stderr, "%s typeinfo differs across DSOs (%p vs %p)\n", name, here, hidden);
  return 1;
}

static int check_throw(void (*thrower)(), const char *what_name) {
  try {
    thrower();
  } catch (const gen::mic::RangeError &err) {
    if (std::strcmp(what_name, "gen::mic::RangeError") == 0 &&
        std::strcmp(err.what(), "cross-dso-range") == 0) {
      return 0;
    }
  } catch (const linkcell::Error &err) {
    if (std::strcmp(what_name, "linkcell::Error") == 0 &&
        std::strcmp(err.what(), "cross-dso-linkcell") == 0) {
      return 0;
    }
  } catch (const pore::BoxError &err) {
    if (std::strcmp(what_name, "pore::BoxError") == 0 &&
        std::strcmp(err.what(), "cross-dso-box") == 0) {
      return 0;
    }
  } catch (const std::exception &err) {
    std::fprintf(stderr, "%s caught std::exception: %s\n", what_name, err.what());
    return 1;
  } catch (...) {
    std::fprintf(stderr, "%s caught an unknown exception\n", what_name);
    return 1;
  }
  std::fprintf(stderr, "%s was not caught by type\n", what_name);
  return 1;
}

int main() {
  int failed = 0;
#if !defined(_WIN32)
  failed |= same("gen::mic::RangeError", &typeid(gen::mic::RangeError), seams_hidden_range_type());
  failed |= same("linkcell::Error", &typeid(linkcell::Error), seams_hidden_linkcell_type());
  failed |= same("pore::BoxError", &typeid(pore::BoxError), seams_hidden_box_type());
#endif
  failed |= check_throw(seams_hidden_throw_range, "gen::mic::RangeError");
  failed |= check_throw(seams_hidden_throw_linkcell, "linkcell::Error");
  failed |= check_throw(seams_hidden_throw_box, "pore::BoxError");
  return failed == 0 ? 0 : 1;
}
