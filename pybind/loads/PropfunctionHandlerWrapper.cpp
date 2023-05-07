// Copyright 2021 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause

#include "PropfunctionHandlerWrapper.h"
#include "pointercollection/pointercollection.h"
#include "pybind11/functional.h"
#include "pybind11/stl.h"

void HierAMuS::PropfunctionHandlerWrapper::registerFunctions() {
  this->temp.def(py::init<>())
      .def("addFunction", &PropfunctionHandler::addFunction)
      .def("getTime", &PropfunctionHandler::getTime)
      .def("set_dt",&PropfunctionHandler::set_dt)
      .def("print",&PropfunctionHandler::print);
}
