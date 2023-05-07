// Copyright 2021 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause

#pragma once
#include <memory>

#include <pybind11/pybind11.h>

namespace py = pybind11;

#include "geometry/GeometryTypes.h"



void GeometryTypesToPybind(pybind11::module &m) {
  py::enum_<HierAMuS::Geometry::GeometryTypes>(m, "GeometryTypes")
      .value("BeamInterface2D",HierAMuS::Geometry::GeometryTypes::BeamInterface2D)
      .value("BeamInterface3D",HierAMuS::Geometry::GeometryTypes::BeamInterface3D)
      .value("Edges",HierAMuS::Geometry::GeometryTypes::Edges)
      .value("Faces",HierAMuS::Geometry::GeometryTypes::Faces)
      .value("Generic",HierAMuS::Geometry::GeometryTypes::Generic)
      .value("LinearEdge",HierAMuS::Geometry::GeometryTypes::LinearEdge)
      .value("LinearBrick",HierAMuS::Geometry::GeometryTypes::LinearBrick)
      .value("LinearTriangle",HierAMuS::Geometry::GeometryTypes::LinearTriangle)
      .value("LinearQuadrilateral",HierAMuS::Geometry::GeometryTypes::LinearQuadrilateral)
      .value("Vertex",
             HierAMuS::Geometry::GeometryTypes::Vertex)

	;
}
