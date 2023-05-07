// Copyright 2021 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause

#include "NodetypesWrapper.h"

void HierAMuS::NodeTypesWrapper::registerFunctions() {
  this->temp.value("displacement", HierAMuS::NodeTypes::displacement)
      .value("rotation", HierAMuS::NodeTypes::rotation)
      .value("undef", HierAMuS::NodeTypes::undef);
}
