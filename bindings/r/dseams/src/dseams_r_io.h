#ifndef DSEAMS_R_IO_H_
#define DSEAMS_R_IO_H_

// Console output of the bundled engine. tools/vendor.sh rewrites the
// engine's std::cout, std::cerr, fprintf(stderr, ...) and std::exit to
// these, so that messages reach the R console and an engine failure
// becomes an R error instead of ending the session.

#include <ostream>

namespace dseams_r {

std::ostream &out();
std::ostream &err();
int eprintf(const char *fmt, ...);
[[noreturn]] void fail(int status);

} // namespace dseams_r

#endif
