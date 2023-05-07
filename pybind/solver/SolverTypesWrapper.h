// Copyright 2021 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause

#pragma once
#include <pybind11/pybind11.h>

namespace py = pybind11;

#include "solver/SolverTypes.h"

namespace HierAMuS {
class SolverTypesWrapper {
public:
  SolverTypesWrapper(py::module &m) : temp(m, "SolverTypes"){};
  void registerFunctions();

private:
  typedef py::enum_<HierAMuS::SolverTypes>
      pw;
  pw temp;
};
} // namespace HierAMuS