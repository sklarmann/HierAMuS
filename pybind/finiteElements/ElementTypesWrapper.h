// Copyright 2021 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause

#pragma once
#include <pybind11/pybind11.h>

namespace py = pybind11;

#include "finiteElements/ElementTypes.h"

namespace HierAMuS {
namespace FiniteElement {

class ElementtypesWrapper {
public:
  ElementtypesWrapper(py::module &m) : temp(m, "Elementtypes"){};
  void registerFunctions();

private:
  typedef py::enum_<Elementtypes> pw;
  pw temp;
};
} // namespace Geometry
} // namespace HierAMuS