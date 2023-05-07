// Copyright 2021 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause

#pragma once
#include <pybind11/pybind11.h>

#include "IntegrationPointsWrapper.h"

namespace py = pybind11;

namespace HierAMuS {
class shapeFunctionsAddClass {
public:
  shapeFunctionsAddClass(py::module &m) : intPoints(m)
           {};
  void registerFunctions();

private:
  IntegrationPointsWrapper intPoints;
};
} // namespace HierAMuS
#undef CNAME