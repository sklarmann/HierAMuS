// Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause


#pragma once

#include <Eigen/SparseCore>
#include <solver/GenericSolver.h>
#ifdef USE_MKL
#include <Eigen/PardisoSupport>
#else
#include <Eigen/OrderingMethods>
#include <Eigen/SparseLU>
#endif

namespace HierAMuS {

class EigenPardisoLLT : public GenericSolver {
public:
  EigenPardisoLLT() = default;
  ;
  ~EigenPardisoLLT() = default;
  ;
  void analyze(Eigen::SparseMatrix<prec, 0, indexType> &SpMat) override;
  void factorize(Eigen::SparseMatrix<prec, 0, indexType> &SpMat) override;
  void solve(Eigen::Matrix<prec, Eigen::Dynamic, 1> &Rhs,
             Eigen::Matrix<prec, Eigen::Dynamic, 1> &solution) override;
  auto solve(Types::MatrixXX<prec> &Rhs) -> Types::MatrixXX<prec> override;

  auto getType() -> SolverTypes override {
    return SolverTypes::TypeEigenPardisoLLT;
  };

private:
#ifdef USE_MKL
  Eigen::PardisoLDLT<Eigen::SparseMatrix<prec, 0, indexType>> solver;
#else
  Eigen::SparseLU<Eigen::SparseMatrix<prec, 0, indexType>,
                  Eigen::AMDOrdering<indexType>>
      solver;
#endif // USE_MKL
};

} // namespace HierAMuS
