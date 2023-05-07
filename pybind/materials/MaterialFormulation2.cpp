// Copyright 2021 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause

#include <memory>

#include <pybind11/pybind11.h>

namespace py = pybind11;

#include "materials/MaterialFormulation2.h"
#include "pointercollection/pointercollection.h"

class PyMaterialFormulation2 : public HierAMuS::Materials::MaterialFormulation2 {
public:
  using HierAMuS::Materials::MaterialFormulation2::MaterialFormulation2;

  //void renew() {
  //  PYBIND11_OVERRIDE(void, HierAMuS::PointerCollection, renew);
  //}
};

void MaterialFormulation2ToPybind(py::module &m) {
  py::class_<HierAMuS::Materials::MaterialFormulation2,
             //PyMaterialFormulation2,
             std::shared_ptr<HierAMuS::Materials::MaterialFormulation2>>(
      m, "MaterialFormulation2")
      .def(py::init<HierAMuS::PointerCollection*>());
}