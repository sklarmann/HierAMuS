// Copyright 2021 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause

#include "EquationHandlerWrapper.h"
#include "pointercollection/pointercollection.h"

namespace HierAMuS {

void EquationHandlerWrapper::registerFunctions() {
  this->temp.def(py::init<PointerCollection*>());
}
} // namespace HierAMuS
