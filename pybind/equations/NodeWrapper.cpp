// Copyright 2021 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause

//
// Created by simon on 1/19/22.
//

#include "NodeWrapper.h"
#include "DegreeOfFreedom.h"
#include <pybind11/stl.h>

namespace HierAMuS {
void GenericNodesWrapper::registerFunctions() {
  this->temp
      .def(
          "getDegreesOfFreedom",
          [](GenericNodes &self) {
            return self.getDegreesOfFreedom();
          },
          py::return_value_policy::reference)
      .def("getDegreeOfFreedom", &GenericNodes::getDegreeOfFreedom,
           py::return_value_policy::reference);
}
} // namespace HierAMuS
