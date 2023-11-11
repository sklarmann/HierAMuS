// Copyright 2021 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause

#pragma once
#include <pybind11/pybind11.h>

namespace py = pybind11;

#include "DegreeOfFreedom.h"

namespace HierAMuS {
//class PyEquationHandler : public EquationHandler {
//public:
//  using GenericFiniteElement::GenericFiniteElement;
//
//  // void renew() {
//  //  PYBIND11_OVERRIDE(void, HierAMuS::PointerCollection, renew);
//  //}
//};

class DegreeOfFreedomWrapper {
public:
  DegreeOfFreedomWrapper(py::module &m) : temp(m, "DegreeOfFreedom"){};
  void registerFunctions();

private:
  typedef py::class_<DegreeOfFreedom, std::shared_ptr<DegreeOfFreedom>>
      pw;
  pw temp;
};
} // namespace HierAMuS