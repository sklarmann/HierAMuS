// Copyright 2021 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <memory>

#include <pybind11/pybind11.h>

namespace py = pybind11;

#include "geometry/Edges.h"

//class PyEdges : public HierAMuS::Geometry::Edges {
//public:
//  using HierAMuS::Geometry::Edges::Edges;
//
//  //void renew() {
//  //  PYBIND11_OVERRIDE(void, HierAMuS::PointerCollection, renew);
//  //}
//};

void GeomEdgesToPybind(py::module &m)
{
  py::class_<HierAMuS::Geometry::Edges,
			 HierAMuS::Geometry::Base,
             std::shared_ptr<HierAMuS::Geometry::Edges>>(m, "Edges")
	;
}