// Copyright 2021 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause

//
// Created by simon on 1/16/22.
//

#include "ScaledBoundary2DWrapper.h"
#include "geometry/ScaledBoundary2D.h"
#include "geometry/GeometryData.h"
#include "pybind11/stl.h"


namespace HierAMuS{
namespace Geometry {

void ScaledBoundary2DWrapper::registerFunctions() {
  this->temp.def("setVerts",&Faces::setVerts)
  .def("setEdges",&Faces::setEdges)
  .def("getVertNums", [](Faces &self){
    std::vector<indexType> vertsOut;
    self.getVerts(vertsOut);
    return vertsOut;
  })
  .def("getEdgeNums",[](Faces &self){
    std::vector<indexType> edgesOut;
    edgesOut.resize(0);
    self.getEdges(edgesOut);
    return edgesOut;
  })
  .def("setScalingCenter",&ScaledBoundary2D::setScalingCenter)
  .def("computeScalingCenter",&ScaledBoundary2D::computeScalingCenter);
}

}
}