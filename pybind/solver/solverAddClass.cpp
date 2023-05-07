// Copyright 2021 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause

#include "solverAddClass.h"

void HierAMuS::solverAdder::registerFunctions()
{
  this->constraints.registerFunctions();
  this->solvers.registerFunctions();
  this->solutions.registerFunctions();

  this->genSolver.registerFunctions();

  this->genSol.registerFunctions();
  this->staticSol.registerFunctions();
  this->staticHomSol.registerFunctions();
  this->newmarkSol.registerFunctions();
}
