// Copyright 2021 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause

#pragma once
#include <memory>

#include <pybind11/pybind11.h>

namespace py = pybind11;

#include "solver/GenericSolutionState.h"
#include "pointercollection/pointercollection.h"

class PyGenericSolutionState: public HierAMuS::GenericSolutionState {
public:
  using HierAMuS::GenericSolutionState::GenericSolutionState;

  //void renew() {
  //  PYBIND11_OVERRIDE(void, HierAMuS::PointerCollection, renew);
  //}
};

void GenericSolutionStateToPybind(pybind11::module &m) {
  py::class_<HierAMuS::GenericSolutionState, PyGenericSolutionState,
             std::shared_ptr<HierAMuS::GenericSolutionState>>(
      m, "GenericSolutionState")
      .def(py::init<HierAMuS::PointerCollection*,HierAMuS::ParameterList&>());
}
