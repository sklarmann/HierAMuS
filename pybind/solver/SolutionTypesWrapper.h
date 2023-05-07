// Copyright 2021 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause

#pragma once
#include <pybind11/pybind11.h>

namespace py = pybind11;

#include "solver/SolutionTypes.h"

namespace HierAMuS {
class SolutionTypesWrapper {
public:
  SolutionTypesWrapper(py::module &m) : temp(m, "SolutionTypes"){};
  void registerFunctions();

private:
  typedef py::enum_<HierAMuS::SolutionTypes> pw;
  pw temp;
};
} // namespace HierAMuS