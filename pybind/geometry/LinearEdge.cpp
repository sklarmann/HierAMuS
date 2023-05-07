// Copyright 2021 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause

#include <memory>

#include "geometry/GeometryData.h"
#include <pybind11/pybind11.h>
#include "plot/vtkplotClass.h"
#include <pybind11/stl.h>

namespace py = pybind11;

#include "geometry/LinearEdge.h"

class PyLinearEdge: public HierAMuS::Geometry::LinearEdge {
public:
  using HierAMuS::Geometry::LinearEdge::LinearEdge;

  //void renew() {
  //  PYBIND11_OVERRIDE(void, HierAMuS::PointerCollection, renew);
  //}
};

void GeomLinearEdgeToPybind(py::module &m)
{
  py::class_<HierAMuS::Geometry::LinearEdge,HierAMuS::Geometry::Edges,
			 HierAMuS::Geometry::Base,
		     PyLinearEdge,
             std::shared_ptr<HierAMuS::Geometry::LinearEdge>>(m, "LinearEdge")
      .def(py::init<>())
      .def("setVerts", &HierAMuS::Geometry::LinearEdge::setVerts)
	;
}