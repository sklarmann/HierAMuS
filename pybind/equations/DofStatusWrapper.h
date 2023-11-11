// Copyright 2021 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause

#pragma once
#include <pybind11/pybind11.h>

namespace py = pybind11;

#include "DofStatus.h"

namespace HierAMuS {

class dofStatusWrapper {
public:
  dofStatusWrapper(py::module &m) : temp(m, "dofStatus"){};
  void registerFunctions();

private:
  typedef py::enum_<dofStatus> pw;
  pw temp;
};
} // namespace HierAMuS