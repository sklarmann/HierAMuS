// Copyright 2021 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause

#include <memory>

#include <pybind11/pybind11.h>

namespace py = pybind11;

#include "geometry/Base.h"

class PyVertex : public HierAMuS::Geometry::Base {
public:
  using HierAMuS::Geometry::Base::Base;

  //void renew() {
  //  PYBIND11_OVERRIDE(void, HierAMuS::PointerCollection, renew);
  //}
};

void BaseToPybind(py::module &m)
{
  py::class_<HierAMuS::Geometry::Base,
             std::shared_ptr<HierAMuS::Geometry::Base>>(m, "Base")
	;
}