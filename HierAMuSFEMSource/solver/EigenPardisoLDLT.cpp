// Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause



#include <iostream>
#include <solver/EigenPardisoLDLT.h>

namespace HierAMuS {

void EigenPardisoLDLT::analyze(Eigen::SparseMatrix<prec, 0, indexType> &SpMat) {
  this->solver.analyzePattern(SpMat);
}

void EigenPardisoLDLT::factorize(
    Eigen::SparseMatrix<prec, 0, indexType> &SpMat) {
  this->solver.factorize(SpMat);
}

void EigenPardisoLDLT::solve(Eigen::Matrix<prec, Eigen::Dynamic, 1> &Rhs,
                             Eigen::Matrix<prec, Eigen::Dynamic, 1> &solution) {
  if (this->solver.info() == Eigen::ComputationInfo::Success) {
    solution = this->solver.solve(Rhs);
  } else {
    solution.setZero();
  }
}
auto EigenPardisoLDLT::solve(Types::MatrixXX<prec> &Rhs)
    -> Types::MatrixXX<prec>  {
  return this->solver.solve(Rhs);
};

} // namespace HierAMuS
