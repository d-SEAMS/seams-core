#include "dseams_r_io.h"

#include <cstdarg>
#include <cstdio>
#include <stdexcept>
#include <streambuf>
#include <string>
#include <thread>
#include <vector>

#include <R_ext/Print.h>

namespace dseams_r {
namespace {

// The library is loaded by dyn.load on R's main thread. R's printing API
// must not be called from any other thread, so text written by an engine
// worker thread is dropped.
const std::thread::id main_thread = std::this_thread::get_id();

bool on_main_thread() { return std::this_thread::get_id() == main_thread; }

void emit(const char *text, std::size_t n, bool to_err) {
  if (n == 0 || !on_main_thread()) {
    return;
  }
  const int len = static_cast<int>(n);
  if (to_err) {
    REprintf("%.*s", len, text);
  } else {
    Rprintf("%.*s", len, text);
  }
}

class ConsoleBuf : public std::streambuf {
public:
  explicit ConsoleBuf(bool to_err) : to_err_(to_err) {}

protected:
  int_type overflow(int_type ch) override {
    if (traits_type::eq_int_type(ch, traits_type::eof())) {
      return traits_type::not_eof(ch);
    }
    const char c = traits_type::to_char_type(ch);
    emit(&c, 1, to_err_);
    return ch;
  }

  std::streamsize xsputn(const char *s, std::streamsize n) override {
    emit(s, static_cast<std::size_t>(n), to_err_);
    return n;
  }

private:
  bool to_err_;
};

} // namespace

std::ostream &out() {
  thread_local ConsoleBuf buf(false);
  thread_local std::ostream stream(&buf);
  return stream;
}

std::ostream &err() {
  thread_local ConsoleBuf buf(true);
  thread_local std::ostream stream(&buf);
  return stream;
}

int eprintf(const char *fmt, ...) {
  va_list args;
  va_start(args, fmt);
  va_list again;
  va_copy(again, args);
  const int n = std::vsnprintf(nullptr, 0, fmt, args);
  va_end(args);
  if (n > 0) {
    std::vector<char> text(static_cast<std::size_t>(n) + 1);
    std::vsnprintf(text.data(), text.size(), fmt, again);
    emit(text.data(), static_cast<std::size_t>(n), true);
  }
  va_end(again);
  return n;
}

void fail(int status) {
  throw std::runtime_error("d-SEAMS engine stopped with status " +
                           std::to_string(status));
}

} // namespace dseams_r
