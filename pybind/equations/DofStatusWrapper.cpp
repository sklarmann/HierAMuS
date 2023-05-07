// Copyright 2021 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause

#include "DofStatusWrapper.h"

void HierAMuS::dofStatusWrapper::registerFunctions() {
  this->temp.value("active", HierAMuS::dofStatus::active)
      .value("inactive", HierAMuS::dofStatus::inactive)
      .value("constraint", HierAMuS::dofStatus::constraint);
}
