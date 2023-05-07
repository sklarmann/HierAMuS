// Copyright 2021 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause

/*
 * GenericSolver.h
 *
 *  Created on: 02.10.2016
 *      Author: Klarmann
 */



#include <Eigen/OrderingMethods>
#include <Eigen/SparseCholesky>
#include <Eigen/SparseCore>
#include <solver/EigenSimplicialLLT.h>

namespace HierAMuS {

void EigenSimplicialLLT::analyze(
    Eigen::SparseMatrix<prec, 0, indexType> &SpMat) {
  this->solver.analyzePattern(SpMat);
}

void EigenSimplicialLLT::factorize(
    Eigen::SparseMatrix<prec, 0, indexType> &SpMat) {
  this->solver.factorize(SpMat);
  // std::cout << "Determinant:  " << this->solver.determinant();
  // std::cout << "Solver Info:  " << this->solver.info() << "   " <<
  // Eigen::Success;
}

void EigenSimplicialLLT::solve(
    Eigen::Matrix<prec, Eigen::Dynamic, 1> &Rhs,
    Eigen::Matrix<prec, Eigen::Dynamic, 1> &solution) {
  if (this->solver.info() == Eigen::ComputationInfo::Success) {
    solution = this->solver.solve(Rhs);
  } else {
    solution.setZero();
  }
};

} // namespace HierAMuS
