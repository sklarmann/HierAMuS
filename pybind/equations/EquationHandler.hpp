// Copyright 2021 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause

#pragma once
#include <memory>

#include <pybind11/pybind11.h>



namespace py = pybind11;

#include "equations/EquationHandler.h"

//class PyEquationHandler : public HierAMuS::EquationHandler {
//public:
//  using HierAMuS::EquationHandler::EquationHandler;
//
//  //void renew() {
//  //  PYBIND11_OVERRIDE(void, HierAMuS::PointerCollection, renew);
//  //}
//};

void EquationHandlerToPybind(pybind11::module &m) {
  py::class_<HierAMuS::EquationHandler,
             std::shared_ptr<HierAMuS::EquationHandler>>(
      m, "EquationHandler")
      .def(py::init<HierAMuS::PointerCollection*>())
      .def("requestNodeSetSetup",&HierAMuS::EquationHandler::requestNodeSetSetup)
      .def("getNumberOfNodeSets",&HierAMuS::EquationHandler::getNumberOfNodeSets)
	;
}
