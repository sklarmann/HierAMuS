// Copyright 2021 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause

#include "VertexWrapper.h"

void HierAMuS::Geometry::VertexWrapper::registerFunctions() {
  this->temp.def(py::init<>())
      .def("setCoordinates",
           py::overload_cast<prec, prec, prec>(
               &HierAMuS::Geometry::Vertex::setCoordinates))
      .def("connectEdge",&HierAMuS::Geometry::Vertex::connectEdge)
      .def("connectFace",&HierAMuS::Geometry::Vertex::connectFace);
}
