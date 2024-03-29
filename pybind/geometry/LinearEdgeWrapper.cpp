// Copyright 2021 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause

#include "geometry/GeometryData.h"

#include "LinearEdgeWrapper.h"
#include "pybind11/stl.h"
#include "plot/vtkplotClass.h"

#include <pointercollection/pointercollection.h>
#include "geometry/GeometryData.h"
#include <geometry/GeometryBaseData.h>

void HierAMuS::Geometry::LinearEdgeWrapper::registerFunctions() {
  this->temp.def(py::init<>())
      .def("setVerts", &HierAMuS::Geometry::LinearEdgeData::setVerts)
      .def("getVerts", &HierAMuS::Geometry::LinearEdgeData::getVertexNumber);
	;
}
