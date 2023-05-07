// Copyright 2021 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause

#pragma once
#include <pybind11/pybind11.h>

namespace py = pybind11;

#include "EdgesWrapper.h"
#include "LinearEdgeWrapper.h"
#include "GeometryDataWrapper.h"
#include "GeometryTypesWrapper.h"
#include "VertexWrapper.h"
#include "BaseWrapper.h"
#include "FacesWrapper.h"
#include "VolumesWrapper.h"
#include "SpecialWrapper.h"
#include "ScaledBoundary2DWrapper.h"



namespace HierAMuS {
class geometryAddClass {
public:
  geometryAddClass(py::module &m)
      : base(m), geoFaces(m), geoVolumes(m), geoSpecial(m), geoData(m), geoTypes(m), geoEdges(m), geoLinEdge(m), vert(m), geoScaledBoundary2D(m)
      {};
  void registerFunctions();

private:
  Geometry::BaseWrapper base;
  Geometry::EdgesWrapper geoEdges;
  Geometry::FacesWrapper geoFaces;
  Geometry::VolumesWrapper geoVolumes;
  Geometry::SpecialWrapper geoSpecial;
  Geometry::VertexWrapper vert;
  Geometry::GeometryDataWrapper geoData;
  Geometry::GeometryTypesWrapper geoTypes;
  Geometry::LinearEdgeWrapper geoLinEdge;
  Geometry::ScaledBoundary2DWrapper geoScaledBoundary2D;
};
} // namespace HierAMuS
#undef CNAME