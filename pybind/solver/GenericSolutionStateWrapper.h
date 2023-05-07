// Copyright 2021 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause

#pragma once
#include <pybind11/pybind11.h>

namespace py = pybind11;

#include "pointercollection/pointercollection.h"
#include "solver/GenericSolutionState.h"

namespace HierAMuS {
class PyGenericSolutionState : public HierAMuS::GenericSolutionState {
public:
  using HierAMuS::GenericSolutionState::GenericSolutionState;

  // void renew() {
  //  PYBIND11_OVERRIDE(void, HierAMuS::PointerCollection, renew);
  //}
};

class GenericSolutionStateWrapper {
public:
  GenericSolutionStateWrapper(py::module &m)
      : temp(m, "GenericSolutionState"){};
  void registerFunctions();

private:
  typedef py::class_<GenericSolutionState, PyGenericSolutionState,
                     std::shared_ptr<GenericSolutionState>>
      pw;
  pw temp;
};
} // namespace HierAMuS