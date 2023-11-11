// Copyright 2021 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause

#include <memory>

#include <pybind11/pybind11.h>

namespace py = pybind11;

#include "geometry/GeometryBaseData.h"

class PyVertex : public HierAMuS::Geometry::GeometryBaseData {
public:
  using HierAMuS::Geometry::GeometryBaseData::GeometryBaseData;

  //void renew() {
  //  PYBIND11_OVERRIDE(void, HierAMuS::PointerCollection, renew);
  //}
};

void BaseToPybind(py::module &m)
{
  py::class_<HierAMuS::Geometry::GeometryBaseData,
             std::shared_ptr<HierAMuS::Geometry::GeometryBaseData>>(m, "GeometryBaseData")
	;
}