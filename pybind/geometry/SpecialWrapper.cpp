// Copyright 2021 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause

//
// Created by simon on 1/16/22.
//

#include "SpecialWrapper.h"
#include "geometry/Special/Special.h"
#include "geometry/GeometryData.h"
#include "pybind11/stl.h"


namespace HierAMuS{
namespace Geometry {

void SpecialWrapper::registerFunctions() {
  this->temp.def("setVerts",&Special::setVerts)
      .def("setEdges",&Special::setEdges)
      .def("setFaces",&Special::setFaces)
      .def("setBeamVertex",&Special::setBeamVertex);
}

}
}