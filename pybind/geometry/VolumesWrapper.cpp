// Copyright 2021 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause

//
// Created by simon on 1/16/22.
//

#include "VolumesWrapper.h"
#include "geometry/Volumes.h"
#include "geometry/GeometryData.h"
#include "pybind11/stl.h"


namespace HierAMuS{
namespace Geometry {

void VolumesWrapper::registerFunctions() {
  this->temp.def("setVerts",&Volumes::setVerts)
      .def("setEdges",&Volumes::setEdges)
      .def("setFaces",&Volumes::setFaces);
}

}
}