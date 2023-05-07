// Copyright 2021 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause

#pragma once
#include <memory>

#include <pybind11/pybind11.h>

namespace py = pybind11;

#include "equations/GenericNodes.h"

// class PyEquationHandler : public HierAMuS::EquationHandler {
// public:
//  using HierAMuS::EquationHandler::EquationHandler;
//
//  //void renew() {
//  //  PYBIND11_OVERRIDE(void, HierAMuS::PointerCollection, renew);
//  //}
//};

void NodeSetToPybind(pybind11::module &m) {
  py::class_<HierAMuS::GenericNodes, std::shared_ptr<HierAMuS::GenericNodes>>(
      m, "GenericNodes")
      .def(py::init<>())
	;
}
