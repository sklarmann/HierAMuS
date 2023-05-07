// Copyright 2021 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause

#pragma once
#include <pybind11/pybind11.h>

namespace py = pybind11;

#include "GenericElementFormulationWrapper.h"
#include "EL209_Python2DWrapper.h"


namespace HierAMuS {
class elementFormulationsAddClass {
public:
  elementFormulationsAddClass(py::module &m) : genElemForm(m), py2DElem(m) {};
  void registerFunctions();

private:
  Elementformulations::GenericElementFormulationWrapper genElemForm;
  Elementformulations::EL290_2DPythonElementWrapper py2DElem;
};
} // namespace HierAMuS