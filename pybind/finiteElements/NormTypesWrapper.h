// Copyright 2021 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause

#pragma once
#include <pybind11/pybind11.h>

namespace py = pybind11;

#include "finiteElements/NormTypes.h"

namespace HierAMuS {
namespace FiniteElement {

class NormtypesWrapper {
public:
  NormtypesWrapper(py::module &m) : temp(m, "NormTypes"){};
  void registerFunctions();

private:
  typedef py::enum_<NormTypes> pw;
  pw temp;
};
} // namespace Geometry
} // namespace HierAMuS