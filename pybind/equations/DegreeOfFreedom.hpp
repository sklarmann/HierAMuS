// Copyright 2021 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause

#pragma once
#include <memory>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;

#include "equations/DegreeOfFreedom.h"

// class PyEquationHandler : public HierAMuS::EquationHandler {
// public:
//  using HierAMuS::EquationHandler::EquationHandler;
//
//  //void renew() {
//  //  PYBIND11_OVERRIDE(void, HierAMuS::PointerCollection, renew);
//  //}
//};

void DegreeOfFreedomToPybind(pybind11::module &m) {
  py::class_<HierAMuS::DegreeOfFreedom,
             std::shared_ptr<HierAMuS::DegreeOfFreedom>>(
      m, "DegreeOfFreedom")
      .def(py::init<>())
	;
}
