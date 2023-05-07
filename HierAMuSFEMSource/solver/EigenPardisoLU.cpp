// Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause


#include <Eigen/SparseCore>
#include <iostream>

#include <solver/EigenPardisoLU.h>
#include <solver/GenericSolver.h>

namespace HierAMuS {

EigenPardisoLU::~EigenPardisoLU() {
  this->first = true;
}

void EigenPardisoLU::analyze(Eigen::SparseMatrix<prec, 0, indexType> &SpMat) {
  
}

void EigenPardisoLU::factorize(Eigen::SparseMatrix<prec, 0, indexType> &SpMat) {
  this->solver.compute(SpMat);
}

void EigenPardisoLU::solve(Eigen::Matrix<prec, Eigen::Dynamic, 1> &Rhs,
                           Eigen::Matrix<prec, Eigen::Dynamic, 1> &solution) {

  if (this->solver.info() == Eigen::ComputationInfo::Success) {
    solution = this->solver.solve(Rhs);
  } else {
    solution.setZero();
  }
}
auto EigenPardisoLU::solve(Types::MatrixXX<prec> &Rhs)
    -> Types::MatrixXX<prec> {
  return this->solver.solve(Rhs);
};

} // namespace HierAMuS
