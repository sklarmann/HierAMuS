// Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause

#pragma once
#include <datatypes.h>

#include <Eigen/Dense>
#include <Eigen/SparseCore>

#include <Base/FEMBase.h>
#include "types/MatrixTypes.h"
#include "SolverTypes.h"

namespace HierAMuS {

class GenericSolver : public FEMBase {
public:
  GenericSolver();
  virtual ~GenericSolver();
  virtual void analyze(Eigen::SparseMatrix<prec, 0, indexType> &SpMat);
  virtual void factorize(Eigen::SparseMatrix<prec, 0, indexType> &SpMat);
  virtual void solve(Eigen::Matrix<prec, Eigen::Dynamic, 1> &Rhs,
                     Eigen::Matrix<prec, Eigen::Dynamic, 1> &solution);
  virtual auto solve(Types::MatrixXX<prec> &Rhs) -> Types::MatrixXX<prec>;

  virtual auto getType() -> SolverTypes { return SolverTypes::Generic; };
};

} // namespace HierAMuS