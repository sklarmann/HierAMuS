// Copyright 2021 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause

#include "controlAddClass.h"

void HierAMuS::controlAddClass::registerFunctions() {
  this->infoData.registerFunctions();
  this->outputData.registerFunctions();
  this->paramList.registerFunctions();
  
}
