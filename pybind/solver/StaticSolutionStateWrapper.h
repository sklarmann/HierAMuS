// Copyright 2021 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause

#pragma once
#include <pybind11/pybind11.h>

namespace py = pybind11;

#include "solver/StaticSolutionState.h"
#include "pointercollection/pointercollection.h"



namespace HierAMuS {
class PyStaticSolutionState : public HierAMuS::StaticSolutionState {
public:
  using HierAMuS::StaticSolutionState::StaticSolutionState;

  // void renew() {
  //  PYBIND11_OVERRIDE(void, HierAMuS::PointerCollection, renew);
  //}
};

class StaticSolutionStateWrapper {
public:
  StaticSolutionStateWrapper(py::module &m)
      : temp(m, "StaticSolutionState"){};
  void registerFunctions();

private:
  typedef py::class_<StaticSolutionState, PyStaticSolutionState,
				     GenericSolutionState,
                     std::shared_ptr<StaticSolutionState>>
      pw;
  pw temp;
};
} // namespace HierAMuS