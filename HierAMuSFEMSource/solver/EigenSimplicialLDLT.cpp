// Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause



#include <Eigen/OrderingMethods>
#include <Eigen/SparseCholesky>
#include <Eigen/SparseCore>
#include <control/HandlingStructs.h>

#include <solver/EigenSimplicialLDLT.h>
#include <solver/GenericSolver.h>

#include <iostream>

namespace HierAMuS {

void EigenSimplicialLDLT::analyze(
    Eigen::SparseMatrix<prec, 0, indexType> &SpMat) {
  this->solver.analyzePattern(SpMat);
}

void EigenSimplicialLDLT::factorize(
    Eigen::SparseMatrix<prec, 0, indexType> &SpMat) {
  this->solver.factorize(SpMat);
  //prec mn = 0, mx = 0;
  //mn = this->solver.vectorD().minCoeff();
  //mx = this->solver.vectorD().maxCoeff();
  // std::cout << "\nConditioning: " << mn << " max " << mx << " rel: " << mx/mn
  // << std::endl; std::cout << "-------------------------------------------" <<
  // std::endl; std::cout << this->solver.vectorD() << std::endl; std::cout <<
  // "-------------------------------------------" << std::endl;
}

void EigenSimplicialLDLT::solve(
    Eigen::Matrix<prec, Eigen::Dynamic, 1> &Rhs,
    Eigen::Matrix<prec, Eigen::Dynamic, 1> &solution) {
  if (this->solver.info() == Eigen::ComputationInfo::Success) {
    solution = this->solver.solve(Rhs);
  } else {
    solution.setZero();
  }
}
auto EigenSimplicialLDLT::solve(Types::MatrixXX<prec> &Rhs)
    -> Types::MatrixXX<prec> {
  return this->solver.solve(Rhs);
};

} // namespace HierAMuS
