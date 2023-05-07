// Copyright 2021 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause

#include "elementFormulationsAddClass.h"

void HierAMuS::elementFormulationsAddClass::registerFunctions()
{
 this->genElemForm.registerFunctions();
 this->py2DElem.registerFunctions();
}
