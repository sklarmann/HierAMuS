// Copyright 2021 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause

#include "GeometryTypesWrapper.h"
#include "geometry/GeometryTypes.h"


void HierAMuS::Geometry::GeometryTypesWrapper::registerFunctions() {
  this->temp.value("Generic", HierAMuS::Geometry::GeometryTypes::Generic)
      .value("Vertex", HierAMuS::Geometry::GeometryTypes::Vertex)

      .value("Edges", HierAMuS::Geometry::GeometryTypes::Edges)
      .value("LinearEdge", HierAMuS::Geometry::GeometryTypes::LinearEdge)
      .value("QuadraticEdge",
             HierAMuS::Geometry::GeometryTypes::QuadraticEdge)

      .value("Faces", HierAMuS::Geometry::GeometryTypes::Faces)
      .value("LinearTriangle",
             HierAMuS::Geometry::GeometryTypes::LinearTriangle)
      .value("LinearQuadrilateral",
             HierAMuS::Geometry::GeometryTypes::LinearQuadrilateral)
      .value("QuadraticQuadrilateral",
             HierAMuS::Geometry::GeometryTypes::QuadraticQuadrilateral)
      .value("ScaledBoundary2D",
             HierAMuS::Geometry::GeometryTypes::ScaledBoundary2D)

      .value("Volumes", HierAMuS::Geometry::GeometryTypes::Volumes)
      .value("LinearBrick", HierAMuS::Geometry::GeometryTypes::LinearBrick)
      .value("LinearPrism", HierAMuS::Geometry::GeometryTypes::LinearPrism)

      .value("BeamInterface2D",
             HierAMuS::Geometry::GeometryTypes::BeamInterface2D)
      .value("BeamInterface3D",
             HierAMuS::Geometry::GeometryTypes::BeamInterface3D);

  this->shapeTypes.value("H0", ShapeFunctionTypes::H0)
      .value("H1", ShapeFunctionTypes::H1)
      .value("HDiv", ShapeFunctionTypes::HDiv)
      .value("HCurl", ShapeFunctionTypes::HCurl);

}


