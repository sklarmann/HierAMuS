// Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause


#include <iostream>
#include <solver/EigenPardisoLLT.h>

namespace HierAMuS {

void EigenPardisoLLT::analyze(Eigen::SparseMatrix<prec, 0, indexType> &SpMat) {
  this->solver.analyzePattern(SpMat);
}

void EigenPardisoLLT::factorize(
    Eigen::SparseMatrix<prec, 0, indexType> &SpMat) {
  this->solver.factorize(SpMat);
  
}

void EigenPardisoLLT::solve(Eigen::Matrix<prec, Eigen::Dynamic, 1> &Rhs,
                            Eigen::Matrix<prec, Eigen::Dynamic, 1> &solution) {
  if (this->solver.info() == Eigen::ComputationInfo::Success) {
    solution = this->solver.solve(Rhs);
  } else {
    solution.setZero();
  }
}
auto EigenPardisoLLT::solve(Types::MatrixXX<prec> &Rhs)
    -> Types::MatrixXX<prec> {
  return this->solver.solve(Rhs);
};

} // namespace HierAMuS
