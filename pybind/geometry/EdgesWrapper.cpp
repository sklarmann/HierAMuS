// Copyright 2021 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause

#include "EdgesWrapper.h"
#include "geometry/Edges.h"
#include "geometry/GeometryData.h"
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>


namespace py = pybind11;

void HierAMuS::Geometry::EdgesWrapper::registerFunctions() {
  //this->temp;
  this->temp//.def(py::init<>())
      .def("setVerts", &HierAMuS::Geometry::Edges::setVerts)
      //.def("getVerts", py::overload_cast<std::vector<indexType>&>(&HierAMuS::Geometry::LinearEdge::getVerts))
      //.def("getVerts",py::overload_cast<HierAMuS::PointerCollection&,std::vector<HierAMuS::Geometry::Base*>&>(&HierAMuS::Geometry::LinearEdge::getVerts))

      .def("getVerts", py::overload_cast<std::vector<indexType> &>(
                           &HierAMuS::Geometry::Edges::getVerts));
	
}
