// Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <map>
#include <solver/GenericSolutionState.h>

#include <vector>

namespace HierAMuS {

class PointerCollection;
class ParameterList;

class TransientSolution : public GenericSolutionState {
public:
  TransientSolution(ParameterList &parameter)
      : GenericSolutionState(parameter){};
  ~TransientSolution(){};
  auto getType() -> SolutionTypes override { return SolutionTypes::Transient; }

private:
};
} /* namespace HierAMuS */
