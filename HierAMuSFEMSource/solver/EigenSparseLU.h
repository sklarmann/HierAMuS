// Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause


#pragma once

#include <Eigen/OrderingMethods>
#include <Eigen/SparseCore>
#include <Eigen/SparseLU>
#include <solver/GenericSolver.h>

namespace HierAMuS {

class EigenSparseLU : public GenericSolver {
public:
  EigenSparseLU() = default;
  ;
  ~EigenSparseLU() override = default;
  ;
  void analyze(Eigen::SparseMatrix<prec, 0, indexType> &SpMat) override;
  void factorize(Eigen::SparseMatrix<prec, 0, indexType> &SpMat) override;
  void solve(Eigen::Matrix<prec, Eigen::Dynamic, 1> &Rhs,
             Eigen::Matrix<prec, Eigen::Dynamic, 1> &solution) override;
  auto solve(Types::MatrixXX<prec> &Rhs) -> Types::MatrixXX<prec> override;

  auto getType() -> SolverTypes override {
    return SolverTypes::TypeEigenSparseLU;
  };

private:
  Eigen::SparseLU<Eigen::SparseMatrix<prec, 0, indexType>,
                  Eigen::COLAMDOrdering<indexType>>
      solver;
};

} // namespace HierAMuS