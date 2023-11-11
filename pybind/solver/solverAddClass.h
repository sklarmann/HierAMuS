// Copyright 2021 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause

#pragma once
#include <pybind11/pybind11.h>

namespace py = pybind11;

#include "GenericSolutionStateWrapper.h"
#include "StaticSolutionStateHomogenizationWrapper.h"
#include "StaticSolutionStateWrapper.h"
#include "TransientSolutionNewmarkWrapper.h"

#include "GenericSolverWrapper.h"

#include "ConstraintWrapper.h"
#include "SolutionTypesWrapper.h"
#include "SolverTypesWrapper.h"

#define CNAME solverAdder

namespace HierAMuS {
class CNAME {
public:
  CNAME(py::module &m)
      : genSol(m), staticSol(m), staticHomSol(m), newmarkSol(m), solvers(m),
        solutions(m), genSolver(m), constraints(m){};
  void registerFunctions();

private:
  GenericSolutionStateWrapper genSol;
  StaticSolutionStateWrapper staticSol;
  StaticSolutionStateHomogenizationWrapper staticHomSol;
  TransientSolutionNewmarkWrapper newmarkSol;

  SolverTypesWrapper solvers;
  SolutionTypesWrapper solutions;

  GenericSolverWrapper genSolver;
  ConstraintWrapper constraints;
};
} // namespace HierAMuS
#undef CNAME