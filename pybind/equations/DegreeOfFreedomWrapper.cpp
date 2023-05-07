// Copyright 2021 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause

//
// Created by simon on 1/19/22.
//

#include "DegreeOfFreedomWrapper.h"

namespace HierAMuS {
void DegreeOfFreedomWrapper::registerFunctions() {
  this->temp
      .def("getId", &DegreeOfFreedom::getId)
      .def("setStatus", &DegreeOfFreedom::setStatus)
      .def("getEqId", &DegreeOfFreedom::getEqId);
  ;
}
} // namespace HierAMuS
