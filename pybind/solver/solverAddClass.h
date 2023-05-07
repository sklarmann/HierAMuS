// Copyright 2021 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause

#pragma once
#include <pybind11/pybind11.h>

namespace py = pybind11;

#include "GenericSolutionStateWrapper.h"
#include "StaticSolutionStateWrapper.h"
#include "TransientSolutionNewmarkWrapper.h"
#include "StaticSolutionStateHomogenizationWrapper.h"

#include "GenericSolverWrapper.h"

#include "SolutionTypesWrapper.h"
#include "SolverTypesWrapper.h"
#include "ConstraintWrapper.h"

#define CNAME solverAdder

namespace HierAMuS {
class CNAME {
public:
  CNAME(py::module &m)
      : constraints(m), genSol(m), solvers(m), solutions(m), staticSol(m),
        newmarkSol(m),
        genSolver(m), staticHomSol(m){};
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