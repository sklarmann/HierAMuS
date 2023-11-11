// Copyright 2021 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause

#include "LoadListWrapper.h"
#include "pybind11/stl.h"
#include "pybind11/functional.h"
#include "pointercollection/pointercollection.h"


void HierAMuS::LoadListWrapper::registerFunctions() {
  this->temp.def(py::init<>())
      .def("print", [](LoadList &self, PointerCollection &pointers) {
        self.print(pointers.getSPDLogger());
      } );
}
