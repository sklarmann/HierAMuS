// Copyright 2021 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause

#pragma once
#include <pybind11/pybind11.h>

namespace py = pybind11;

#include "solver/StaticSolutionStateHomogenization.h"
#include "pointercollection/pointercollection.h"



namespace HierAMuS {
class PyStaticSolutionStateHomogenization
    : public HierAMuS::StaticSolutionStateHomogenization {
public:
  using HierAMuS::StaticSolutionStateHomogenization::StaticSolutionStateHomogenization;

  // void renew() {
  //  PYBIND11_OVERRIDE(void, HierAMuS::PointerCollection, renew);
  //}
};

class StaticSolutionStateHomogenizationWrapper {
public:
  StaticSolutionStateHomogenizationWrapper
    (py::module &m)
      : temp(m, "StaticSolutionStateHomogenization"){};
  void registerFunctions();

private:
  typedef py::class_<StaticSolutionStateHomogenization, PyStaticSolutionStateHomogenization,
      StaticSolutionState,
      GenericSolutionState,
                     std::shared_ptr<StaticSolutionStateHomogenization>>
      pw;
  pw temp;
};
} // namespace HierAMuS