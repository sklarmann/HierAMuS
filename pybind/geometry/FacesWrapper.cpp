// Copyright 2021 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause

//
// Created by simon on 1/16/22.
//

#include "FacesWrapper.h"
#include "geometry/Faces/FacesData.h"
#include "geometry/GeometryData.h"
#include "pybind11/stl.h"

namespace HierAMuS {
namespace Geometry {

void FacesWrapper::registerFunctions() {
  this->temp.def("setVerts", &FacesData::setVerts)
      .def("setEdges", &FacesData::setEdges)
      .def("getVertNums", [](FacesData &self) { return self.getVertexNumbers(); })
      .def("getEdgeNums", [](FacesData &self) {
        auto edgesOut = self.getEdgeNumbers();
        return edgesOut;
      });
}

} // namespace Geometry
} // namespace HierAMuS