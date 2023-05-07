// Copyright 2021 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause

//
// Created by simon on 1/19/22.
//

#include "ConstraintWrapper.h"
#include <pybind11/stl.h>

namespace HierAMuS {
void ConstraintWrapper::registerFunctions() {
  this->ConstraintTypesModule.value("GeneralLink", ConstraintTypes::GeneralLink);
  this->GeneralLinkModule.def("set", &GeneralLink::set);

  this->ConstraintHandlerModule.def(
      "GeneralLinkGeo", &ConstraintHandler::GeneralLinkGeo);

}
}
