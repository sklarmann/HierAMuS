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
#include <Eigen/SparseCore>
#include <Eigen/SparseQR>
#include <solver/EigenSparseQR.h>

namespace HierAMuS {

void EigenSparseQR::analyze(Eigen::SparseMatrix<prec, 0, indexType> &SpMat) {
  this->solver.analyzePattern(SpMat);
}

void EigenSparseQR::factorize(Eigen::SparseMatrix<prec, 0, indexType> &SpMat) {
  this->solver.factorize(SpMat);
}

void EigenSparseQR::solve(Eigen::Matrix<prec, Eigen::Dynamic, 1> &Rhs,
                          Eigen::Matrix<prec, Eigen::Dynamic, 1> &solution) {
  if (this->solver.info() == Eigen::ComputationInfo::Success) {
    solution = this->solver.solve(Rhs);
  } else {
    solution.setZero();
  }
};

} // namespace HierAMuS
