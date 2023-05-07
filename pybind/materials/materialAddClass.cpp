// Copyright 2021 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause

#include "materialAddClass.h"


void HierAMuS::materialAddClass::registerFunctions()
{
  this->matlist.registerFunctions();
  this->material.registerFunctions();
  this->elemList.registerFunctions();
  this->matFormList.registerFunctions();
  this->GenericMat.registerFunctions();
}
