#ifndef SEAMS_NUMBER_PARSE_HPP_
#define SEAMS_NUMBER_PARSE_HPP_

#include <cerrno>
#include <cmath>
#include <cstdlib>
#include <locale.h>
#include <string>
#include <string_view>

namespace seams {

// libc++ from the conda macOS 11 toolchain deletes std::from_chars for
// a floating value. strtod_l in the C locale follows the same grammar:
// no leading plus, no hexadecimal form, the whole text consumed, and
// overflow or underflow to zero rejected.

inline bool parseExactDouble(std::string_view text, double &out) {
  if (text.empty() || text.front() == '+') {
    return false;
  }
  for (char raw : text) {
    const auto c = static_cast<unsigned char>(raw);
    if (raw == 'x' || raw == 'X' || c == ' ' || c == '\t' || c == '\n' ||
        c == '\v' || c == '\f' || c == '\r') {
      return false;
    }
  }
#if defined(_WIN32)
  // MinGW and MSVC spell the C-locale converter this way. newlocale
  // is not in their headers. The grammar stays the one strtod_l uses.
  static const _locale_t loc = _create_locale(LC_ALL, "C");
  if (loc == nullptr) {
    return false;
  }
#else
  static const locale_t loc =
      newlocale(LC_ALL_MASK, "C", static_cast<locale_t>(nullptr));
  if (loc == static_cast<locale_t>(nullptr)) {
    return false;
  }
#endif
  const std::string token(text);
  char *end = nullptr;
  errno = 0;
#if defined(_WIN32)
  const double value = _strtod_l(token.c_str(), &end, loc);
#else
  const double value = strtod_l(token.c_str(), &end, loc);
#endif
  if (end != token.c_str() + token.size()) {
    return false;
  }
  if (errno == ERANGE && (value == 0.0 || !std::isfinite(value))) {
    return false;
  }
  out = value;
  return true;
}

} // namespace seams

#endif
