// Copyright 2021 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause

#include "loadsAddClass.h"

void HierAMuS::loadsAddClass::registerFunctions() {
  this->func.registerFunctions();
  this->props.registerFunctions();
  this->loads.registerFunctions();
}
