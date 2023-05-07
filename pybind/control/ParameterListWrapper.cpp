// Copyright 2021 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause

#include "ParameterListWrapper.h"
#include "pybind11/eigen.h"

void HierAMuS::ParameterListWrapper::registerFunctions() {
  this->temp.def(py::init<>())
      .def("add", py::overload_cast<const std::string &, prec>(&ParameterList::add))
      .def("add", py::overload_cast<const std::string &, Types::MatrixXX<prec>>(&ParameterList::add))
	;
}
