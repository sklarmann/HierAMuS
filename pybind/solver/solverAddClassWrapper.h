// Copyright 2021 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause

#pragma once
#include <pybind11/pybind11.h>

namespace py = pybind11;

#include "GenericSolutionStateWrapper.h"

namespace HierAMuS {
class solverAdder {
public:
  solverAdder(py::module &m) : genSol(m){};
  void registerFunctions();

private:
  GenericSolutionStateWrapper genSol;
};
} // namespace HierAMuS