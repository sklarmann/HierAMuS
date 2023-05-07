// Copyright 2021 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause

#pragma once
#include <pybind11/pybind11.h>

namespace py = pybind11;

#include "geometry/GeometryTypes.h"

namespace HierAMuS {
namespace Geometry {

class GeometryTypesWrapper {
public:
  GeometryTypesWrapper(py::module &m)
      : temp(m, "GeometryTypes"), shapeTypes(m, "ShapeFunctionTypes"){};
  void registerFunctions();

private:
  typedef py::enum_<GeometryTypes> pw;
  pw temp;
  py::enum_<ShapeFunctionTypes> shapeTypes;
};
} // namespace Geometry
} // namespace HierAMuS