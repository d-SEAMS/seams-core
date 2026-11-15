// Copyright (c) d-SEAMS developers
// SPDX-License-Identifier: MIT
// Some PLUMED headers declare addComponent(name, shape=empty) while the
// library exports only addComponent(name). Calling the header declaration
// asks the linker for a symbol that is not there. These definitions name
// the one-argument symbols the library does export.

#include "plumed/config/version.h"

#include <string>
#include <vector>

namespace PLMD {

class ActionWithValue {
 public:
#if PLUMED_VERSION_MAJOR > 2 || (PLUMED_VERSION_MAJOR == 2 && PLUMED_VERSION_MINOR >= 10)
  // 2.10 exports the shape overload and not the one-argument symbol.
  void addComponent(const std::string &name, const std::vector<unsigned> &shape);
  void addComponentWithDerivatives(const std::string &name,
                                   const std::vector<unsigned> &shape);
#else
  void addComponent(const std::string &name);
  void addComponentWithDerivatives(const std::string &name);
#endif
};

void addPlainComponent(ActionWithValue *self, const std::string &name) {
#if PLUMED_VERSION_MAJOR > 2 || (PLUMED_VERSION_MAJOR == 2 && PLUMED_VERSION_MINOR >= 10)
  self->addComponent(name, std::vector<unsigned>());
#else
  self->addComponent(name);
#endif
}

void addDerivativeComponent(ActionWithValue *self, const std::string &name) {
#if PLUMED_VERSION_MAJOR > 2 || (PLUMED_VERSION_MAJOR == 2 && PLUMED_VERSION_MINOR >= 10)
  self->addComponentWithDerivatives(name, std::vector<unsigned>());
#else
  self->addComponentWithDerivatives(name);
#endif
}

}  // namespace PLMD
