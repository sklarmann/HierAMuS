// Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause

#include <solver/GenericSolver.h>

namespace HierAMuS {

GenericSolver::GenericSolver() {}

GenericSolver::~GenericSolver() {}

void GenericSolver::analyze(Eigen::SparseMatrix<prec, 0, indexType> &SpMat) {}

void GenericSolver::factorize(Eigen::SparseMatrix<prec, 0, indexType> &SpMat) {}

void GenericSolver::solve(Eigen::Matrix<prec, Eigen::Dynamic, 1> &Rhs,
                          Eigen::Matrix<prec, Eigen::Dynamic, 1> &solution) {}

auto GenericSolver::solve(Types::MatrixXX<prec> &Rhs) -> Types::MatrixXX<prec> {
  return Types::MatrixXX<prec>();
}


} // namespace HierAMuS
