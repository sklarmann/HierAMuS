// Copyright 2021 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause

/*
 * GenericSolver.h
 *
 *  Created on: 02.10.2016
 *      Author: Klarmann
 */

#include <Eigen/SparseCore>
#include <iostream>

#include <solver/EigenPardisoLU.h>
#include <solver/GenericSolver.h>

namespace HierAMuS {

EigenPardisoLU::~EigenPardisoLU() {
  this->first = true;
  // Eigen::SparseMatrix<prec, 0, indexType> SpMat;
  // SpMat.resize(0, 0);
  // this->solver.analyzePattern(SpMat);
  // std::cout << "Pardiso LU destructor called" << std::endl;
}

void EigenPardisoLU::analyze(Eigen::SparseMatrix<prec, 0, indexType> &SpMat) {
  // this->solver.compute(SpMat);
}

void EigenPardisoLU::factorize(Eigen::SparseMatrix<prec, 0, indexType> &SpMat) {
  // if (first) {
  //	this->solver.compute(SpMat);
  //	this->first = false;
  //	std::cout << "frist" << std::endl;
  //}
  // else
  //{
  //	this->solver.factorize(SpMat);
  //	std::cout << "second" << std::endl;
  //}
  this->solver.compute(SpMat);
}

void EigenPardisoLU::solve(Eigen::Matrix<prec, Eigen::Dynamic, 1> &Rhs,
                           Eigen::Matrix<prec, Eigen::Dynamic, 1> &solution) {

  if (this->solver.info() == Eigen::ComputationInfo::Success) {
    solution = this->solver.solve(Rhs);
  } else {
    solution.setZero();
  }
};

} // namespace HierAMuS
