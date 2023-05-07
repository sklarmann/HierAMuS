// Copyright 2021 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause

#include "plotAddClass.h"

void HierAMuS::plotAddClass::registerFunctions() {
  this->vtkplot.registerFunctions();
  this->plotcontrol.registerFunctions();
}
