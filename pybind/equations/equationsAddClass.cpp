// Copyright 2021 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause

#include "equationsAddClass.h"

void HierAMuS::equationsAddClass::registerFunctions() {
  this->nodeTypes.registerFunctions();
  this->eqHandler.registerFunctions();
  this->dofStatus.registerFunctions();
  this->Dofs.registerFunctions();
  this->nodes.registerFunctions();
}
