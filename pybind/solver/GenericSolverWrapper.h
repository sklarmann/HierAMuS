// Copyright 2021 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause

#pragma once
#include <pybind11/pybind11.h>

namespace py = pybind11;

#include "solver/GenericSolver.h"

namespace HierAMuS {
class PyGenericSolver : public HierAMuS::GenericSolver {
public:
  using HierAMuS::GenericSolver::GenericSolver;

  // void renew() {
  //  PYBIND11_OVERRIDE(void, HierAMuS::PointerCollection, renew);
  //}
};

class GenericSolverWrapper {
public:
  GenericSolverWrapper(py::module &m)
      : temp(m, "GenericSolver"){};
  void registerFunctions();

private:
  typedef py::class_<GenericSolver, PyGenericSolver,
                     std::shared_ptr<GenericSolver>>
      pw;
  pw temp;
};
} // namespace HierAMuS