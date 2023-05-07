// Copyright 2021 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause

#include <memory>

#include <pybind11/pybind11.h>

namespace py = pybind11;

#include "materials/MaterialformulationList.h"

class PyMaterialFormulationList : public HierAMuS::Materials::MaterialFormulationList {
public:
  using HierAMuS::Materials::MaterialFormulationList::MaterialFormulationList;

  //void renew() {
  //  PYBIND11_OVERRIDE(void, HierAMuS::PointerCollection, renew);
  //}
};

void MaterialFormulationListToPybind(py::module &m) {
  py::class_<HierAMuS::Materials::MaterialFormulationList,
             //PyMaterialFormulationList,
             std::shared_ptr<HierAMuS::Materials::MaterialFormulationList>>(
      m, "MaterialFormulationList")
      .def(py::init<>());
}