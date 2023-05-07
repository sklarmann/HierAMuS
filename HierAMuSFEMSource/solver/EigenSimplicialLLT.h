// Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause



#pragma once

#include <Eigen/OrderingMethods>
#include <Eigen/SparseCholesky>
#include <Eigen/SparseCore>
#include <solver/GenericSolver.h>

namespace HierAMuS {

class EigenSimplicialLLT : public GenericSolver {
public:
  EigenSimplicialLLT() = default;
  ;
  ~EigenSimplicialLLT() override = default;
  ;
  void analyze(Eigen::SparseMatrix<prec, 0, indexType> &SpMat) override;
  void factorize(Eigen::SparseMatrix<prec, 0, indexType> &SpMat) override;
  void solve(Eigen::Matrix<prec, Eigen::Dynamic, 1> &Rhs,
             Eigen::Matrix<prec, Eigen::Dynamic, 1> &solution) override;
  auto solve(Types::MatrixXX<prec> &Rhs) -> Types::MatrixXX<prec> override;

  auto getType() -> SolverTypes override {
    return SolverTypes::TypeEigenSimplicialLLT;
  };

private:
  Eigen::SimplicialLLT<Eigen::SparseMatrix<prec, 0, indexType>, Eigen::Lower,
                       Eigen::AMDOrdering<indexType>>
      solver;
};

} // namespace HierAMuS
