// Copyright 2021 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause

#include "LinearEdgeWrapper.h"
#include "geometry/GeometryData.h"
#include "pybind11/stl.h"
#include "plot/vtkplotClass.h"

#include <pointercollection/pointercollection.h>
#include <geometry/Base.h>

void HierAMuS::Geometry::LinearEdgeWrapper::registerFunctions() {
  this->temp.def(py::init<>())
      .def("setVerts", &HierAMuS::Geometry::LinearEdge::setVerts)
      //.def("getVerts", py::overload_cast<std::vector<indexType>&>(&HierAMuS::Geometry::LinearEdge::getVerts))
      //.def("getVerts",py::overload_cast<HierAMuS::PointerCollection&,std::vector<HierAMuS::Geometry::Base*>&>(&HierAMuS::Geometry::LinearEdge::getVerts))
      .def("getVerts",
           ([](HierAMuS::Geometry::LinearEdge &self)
      { std::vector<indexType> vec;
             self.getVerts(vec);
             return vec;
           }))
      .def("getVerts", py::overload_cast<std::vector<indexType> &>(
                           &HierAMuS::Geometry::LinearEdge::getVerts));
	;
}
