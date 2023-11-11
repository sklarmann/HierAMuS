// Copyright 2021 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause

#include "VertexWrapper.h"
#include "pybind11/eigen.h"

void HierAMuS::Geometry::VertexWrapper::registerFunctions() {
  this->temp.def(py::init<>())
      .def("setCoordinates",
           py::overload_cast<prec, prec, prec>(
                                 &HierAMuS::Geometry::VertexData::setCoordinates))
      .def("getCoordinates", py::overload_cast<>(&VertexData::getCoordinates))
      .def("connectEdge",&HierAMuS::Geometry::VertexData::connectEdge)
      .def("connectFace",&HierAMuS::Geometry::VertexData::connectFace);
}
