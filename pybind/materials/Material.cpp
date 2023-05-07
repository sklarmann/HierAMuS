// Copyright 2021 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause

#include <memory>

#include <pybind11/pybind11.h>

namespace py = pybind11;

#include "materials/Material.h"
#include "pointercollection/pointercollection.h"

class PyMaterial : public HierAMuS::Materials::Material {
public:
  using HierAMuS::Materials::Material::Material;

  //void renew() {
  //  PYBIND11_OVERRIDE(void, HierAMuS::PointerCollection, renew);
  //}
};

void MaterialToPybind(py::module &m) {
  py::class_<HierAMuS::Materials::Material,
			 //PyMaterial,
             std::shared_ptr<HierAMuS::Materials::Material>>(
      m, "Material")
      .def(py::init<indexType>());
}