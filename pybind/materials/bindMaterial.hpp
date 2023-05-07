// Copyright 2021 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <pybind11/pybind11.h>

namespace py = pybind11;

void ElementFormulationListToPybind(py::module &m);
void GenericMaterialFormulationToPybind(py::module &m);
void MaterialToPybind(py::module &m);
void MaterialFormulation2ToPybind(py::module &m);

void addMaterialFolder(py::module &m) {
  auto MaterialModule = m.def_submodule("Materials", "Material Module");

  ElementFormulationListToPybind(MaterialModule);
  GenericMaterialFormulationToPybind(MaterialModule);
  MaterialToPybind(MaterialModule);
  MaterialFormulation2ToPybind(MaterialModule);
}