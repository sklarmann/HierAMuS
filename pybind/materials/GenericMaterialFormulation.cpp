// Copyright 2021 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause

#include <memory>

#include <pybind11/pybind11.h>

namespace py = pybind11;

#include <materials/GenericMaterialFormulation.h>
#include <pointercollection/pointercollection.h>

class PyGenericMaterialFormulation : public HierAMuS::Materials::GenericMaterialFormulation {
public:
  using HierAMuS::Materials::GenericMaterialFormulation::GenericMaterialFormulation;

  //void renew() {
  //  PYBIND11_OVERRIDE(void, HierAMuS::PointerCollection, renew);
  //}
};

void GenericMaterialFormulationToPybind(py::module &m) {
  py::class_<HierAMuS::Materials::GenericMaterialFormulation,
             PyGenericMaterialFormulation,
      std::shared_ptr<HierAMuS::Materials::GenericMaterialFormulation>>(
      m, "GenericMaterialFormulation")
      .def(py::init<HierAMuS::PointerCollection*>());
}