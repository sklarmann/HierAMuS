// Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <solver/SolutionTypes.h>

#include <Eigen/Sparse>
#include <solver/GenericSolutionState.h>
#include <vector>



namespace HierAMuS {

 class PointerCollection;


class StaticSolutionState : public GenericSolutionState {
public:
  StaticSolutionState(ParameterList &parameter);
  StaticSolutionState(const StaticSolutionState &other);
  virtual ~StaticSolutionState() override;

  auto getCopy() -> std::shared_ptr<GenericSolutionState> override;
  
  auto getType() -> SolutionTypes override { return SolutionTypes::StaticSolutionState; }

  void setInitialValues(indexType numberOfEquations,
                        indexType numberOfActiveEquations) override;
  
  void setSparseMatrix(PointerCollection& pointers) override;
  void insertStiffnessResidual(
      Eigen::Matrix<prec, Eigen::Dynamic, Eigen::Dynamic> &stiffness,
      Eigen::Matrix<prec, Eigen::Dynamic, 1> &residual,
      std::vector<DegreeOfFreedom *> &Dofs) override;

  void computeLoads(PointerCollection &pointers) override;
  void assembleSystem(PointerCollection& pointers) override;

  
  void nextSolutionStep() override;

  void factorize() override;
  void solve(PointerCollection& pointers) override;
  auto lgsResidual() -> prec override;;
  auto energyNorm() -> prec override { return this->energyVal; };
  auto residual() -> prec override { return this->residualVal; };
  void setEquationZero() override;
  void updateSolution(PointerCollection& pointers) override;
  void dampedSolutionUpdate(PointerCollection& pointers) override;
  void computeEigenValues(PointerCollection& pointers, indexType number,
                          indexType addNumber = 0, bool max = false,
                          prec tol = 1e-10, prec shift = 1e-10) override;


  auto getSolution(GenericNodes &node) -> Types::Vector3<prec> override;
  auto getSolution(indexType globalId) -> prec override;
  auto getSolution(std::vector<DegreeOfFreedom*>& Dofs) -> Eigen::Matrix<prec, Eigen::Dynamic, 1> override;
  auto getIncrementalSolution(std::vector<DegreeOfFreedom*>& Dofs) -> Eigen::Matrix<prec, Eigen::Dynamic, 1> override;
  auto getNewtonSolution(std::vector<DegreeOfFreedom*>& Dofs) -> Eigen::Matrix<prec, Eigen::Dynamic, 1> override;


  void resetSolution() override;

  void printSpMat(PointerCollection& pointers) override;

  void computeConditionNumber(PointerCollection& pointers) override;

  void ctestout(PointerCollection& pointers) override;

  
  void toFile(PointerCollection &pointers, std::ofstream &out) override;
  void fromFile(PointerCollection &pointers, std::ifstream &in) override;

  
  void RVEDatatoFile(PointerCollection &pointers, std::ofstream &out) override;
  void RVEDatafromFile(PointerCollection &pointers, std::ifstream &in) override;
  

protected:
  Eigen::Matrix<prec, Eigen::Dynamic, 1> Solution, IncSolution, dIncSolution, NewtonSolution;

  //		Eigen::SparseLU<Eigen::SparseMatrix<prec, 0, indexType>,
  //Eigen::COLAMDOrdering<indexType> > solver;

  Types::VectorX<prec> Rhs, RhsB;
  Types::VectorX<prec> eqSol, incSol;
  std::vector<Eigen::Triplet<prec,indexType>> tripletList;
  Types::SparseMatrix<prec,indexType> Kaa, Kab, Kbb, Kba;

private:
  void insertStiffnessResidualSymUpper(
      Eigen::Matrix<prec, Eigen::Dynamic, Eigen::Dynamic> &stiffness,
      Eigen::Matrix<prec, Eigen::Dynamic, 1> &residual,
      std::vector<DegreeOfFreedom *> &Dofs);
  void insertStiffnessResidualSymLower(
      Eigen::Matrix<prec, Eigen::Dynamic, Eigen::Dynamic> &stiffness,
      Eigen::Matrix<prec, Eigen::Dynamic, 1> &residual,
      std::vector<DegreeOfFreedom *> &Dofs);
  void insertStiffnessResidualFull(
      Eigen::Matrix<prec, Eigen::Dynamic, Eigen::Dynamic> &stiffness,
      Eigen::Matrix<prec, Eigen::Dynamic, 1> &residual,
      std::vector<DegreeOfFreedom *> &Dofs);


  
  void computeEigenValuesSymUpper(PointerCollection &pointers, indexType number,
                                  indexType addNumber, bool max,
                                  prec tol, prec shift);
  void computeEigenValuesSymLower(PointerCollection &pointers, indexType number,
                                  indexType addNumber, bool max,
                                  prec tol, prec shift);
  void computeEigenValuesUnsym(PointerCollection &pointers, indexType number,
                                  indexType addNumber, bool max,
                                  prec tol, prec shift);

};
} /* namespace HierAMuS */
