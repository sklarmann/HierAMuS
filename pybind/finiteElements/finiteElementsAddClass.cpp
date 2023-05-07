// Copyright 2021 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause

#include "finiteElementsAddClass.h"

void HierAMuS::finiteElementsAddClass::registerFunctions()
{
  this->elemList.registerFunctions();
  this->genericElem.registerFunctions();
  this->elemTypes.registerFunctions();
  this->normTypes.registerFunctions();
}
