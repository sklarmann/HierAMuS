// Copyright 2021 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause

#pragma once
#include <pybind11/pybind11.h>

namespace py = pybind11;

#include "equations/EquationHandler.h"

namespace HierAMuS {
//class PyEquationHandler : public EquationHandler {
//public:
//  using GenericFiniteElement::GenericFiniteElement;
//
//  // void renew() {
//  //  PYBIND11_OVERRIDE(void, HierAMuS::PointerCollection, renew);
//  //}
//};

class EquationHandlerWrapper {
public:
  EquationHandlerWrapper(py::module &m) : temp(m, "EquationHandler"){};
  void registerFunctions();

private:
  typedef py::class_<EquationHandler, std::shared_ptr<EquationHandler>>
      pw;
  pw temp;
};
} // namespace HierAMuS