// Copyright 2021 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause

#include "TimefunctionWrapper.h"
#include "pybind11/stl.h"
#include "pybind11/functional.h"

void HierAMuS::TimefunctionWrapper::registerFunctions() {
  this->fun.def(py::init<>())
      .def("setLambda", &HierAMuS::Function::setLambda)
      .def("evaluate", &HierAMuS::Function::evaluate);
}
