// Copyright (c) d-SEAMS developers
// SPDX-License-Identifier: MIT
// Some PLUMED headers declare addComponent(name, shape=empty) while the
// library exports only addComponent(name). Calling the header declaration
// asks the linker for a symbol that is not there. These definitions name
// the one-argument symbols the library does export.

#include <string>

namespace PLMD {

class ActionWithValue {
 public:
  void addComponent(const std::string &name);
  void addComponentWithDerivatives(const std::string &name);
};

void addPlainComponent(ActionWithValue *self, const std::string &name) {
  self->addComponent(name);
}

void addDerivativeComponent(ActionWithValue *self, const std::string &name) {
  self->addComponentWithDerivatives(name);
}

}  // namespace PLMD
