// Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <map>
#include <solver/TransientSolution.h>

#include <vector>

namespace HierAMuS {

class PointerCollection;

class TransientSolutionNewmark : public TransientSolution {
public:
  TransientSolutionNewmark(ParameterList &parameter);
  ~TransientSolutionNewmark() override;
  void setValues(std::map<std::string, prec> &values) override;
  void setInitialValues(indexType numberOfEquations,
                        indexType numberOfActiveEquations) override;

  void setSparseMatrix(PointerCollection &pointers) override;
  void insertStiffnessResidual(
      Eigen::Matrix<prec, Eigen::Dynamic, Eigen::Dynamic> &stiffness,
      Eigen::Matrix<prec, Eigen::Dynamic, 1> &residual,
      std::vector<DegreeOfFreedom *> &Dofs) override;

  void computeLoads(PointerCollection &pointers) override;

  void assembleSystem(PointerCollection &pointers) override;
  void setEquationZero() override;

  Eigen::Matrix<prec, Eigen::Dynamic, 1>
  getSolution(std::vector<DegreeOfFreedom *> &Dofs) override;
  Eigen::Matrix<prec, Eigen::Dynamic, 1>
  getVelocity(std::vector<DegreeOfFreedom *> Dofs) override;
  Eigen::Matrix<prec, Eigen::Dynamic, 1>
  getAcceleration(std::vector<DegreeOfFreedom *> Dofs) override;
  void factorize() override;
  void solve(PointerCollection &pointers) override;
  // void setEquationZero();
  void updateSolution(PointerCollection &pointers) override;
  // void computeEigenValues();
  //
  prec getSolution(indexType globalId) override;

  void nextSolutionStep() override;

private:
  void insertMassMatrix(
      Eigen::Matrix<prec, Eigen::Dynamic, Eigen::Dynamic> &stiffness,
      std::vector<DegreeOfFreedom *> &Dofs);
  prec beta, gamma;

  Eigen::Matrix<prec, Eigen::Dynamic, 1> Solution, IncSolution, dIncSolution;
  Eigen::Matrix<prec, Eigen::Dynamic, 1> vn, vn1, an, an1;
  Eigen::Matrix<prec, Eigen::Dynamic, 1> compW, uncompW;
  Eigen::SparseMatrix<prec, 0, indexType> Stiffness, Mass, Damping;
  Eigen::SparseMatrix<prec, 0, indexType> SpMat;

  Eigen::Matrix<prec, Eigen::Dynamic, 1> Rhs;
  Eigen::Matrix<prec, Eigen::Dynamic, 1> eqSol;
  std::vector<Eigen::Triplet<prec, indexType>> tripletList;

  //		Eigen::SparseLU<Eigen::SparseMatrix<prec, 0, indexType>,
  //Eigen::COLAMDOrdering<indexType> > solver;

  // Eigen::Matrix<prec, Eigen::Dynamic, 1> Rhs;
  // Eigen::Matrix<prec, Eigen::Dynamic, 1> eqSol;
  // std::vector<Eigen::Triplet> tripletList;
};
} /* namespace HierAMuS */
