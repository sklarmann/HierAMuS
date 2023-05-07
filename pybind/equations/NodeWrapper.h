// Copyright 2021 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause

#pragma once
#include <pybind11/pybind11.h>

namespace py = pybind11;

#include "equations/GenericNodes.h"

namespace HierAMuS {
//class PyEquationHandler : public EquationHandler {
//public:
//  using GenericFiniteElement::GenericFiniteElement;
//
//  // void renew() {
//  //  PYBIND11_OVERRIDE(void, HierAMuS::PointerCollection, renew);
//  //}
//};

class GenericNodesWrapper {
public:
  GenericNodesWrapper(py::module &m) : temp(m, "GenericNodes"){};
  void registerFunctions();

private:
  typedef py::class_<GenericNodes, std::shared_ptr<GenericNodes>>
      pw;
  pw temp;
};
} // namespace HierAMuS