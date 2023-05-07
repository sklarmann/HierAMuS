// Copyright 2021 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause

#include "SolutionTypesWrapper.h"

void HierAMuS::SolutionTypesWrapper::registerFunctions() {
  this->temp
      .value("GenericSolutionState",
             HierAMuS::SolutionTypes::GenericSolutionState)
      .value("LinearStaticSolutionState",
             HierAMuS::SolutionTypes::LinearStaticSolutionState)
      .value("StaticSolutionState",
             HierAMuS::SolutionTypes::StaticSolutionState)
      .value("StaticSolutionHomogenization",
             HierAMuS::SolutionTypes::StaticSolutionHomogenization)
      .value("Transient", HierAMuS::SolutionTypes::Transient)
      .value("TransientSolutionNewmark",
             HierAMuS::SolutionTypes::TransientSolutionNewmark);
}
