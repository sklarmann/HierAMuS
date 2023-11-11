// Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause

#pragma once


#include <solver/SolutionTypes.h>
#include <solver/SolverTypes.h>

#include "HistoryDataNew/HistoryDataIterator.h"
#include "HistoryDataNew/HistoryDataManager.h"

#include "Constraints/ConstraintHandler.h"

#include <Eigen/SparseCore>

#include <Eigen/Dense>

#include <types/MatrixTypes.h>

namespace HierAMuS {

class GenericNodes;
class GenericSolver;
enum class HistoryTypes { Constant, TimeUpdate };
class ParameterList;
class PropfunctionHandler;

class ElementDataFields {

  public:
  void request_field(indexType elementId, indexType fieldId, indexType rows,
                     indexType cols);
  auto get_field(indexType elementId, indexType fieldId)
      -> Types::MatrixXX<prec> &;
  void set_field(indexType elementId, indexType fieldId,
                 Types::MatrixXX<prec> &data);

  private:
  std::map<indexType, std::vector<Types::MatrixXX<prec>>> m_data;
};


class GenericSolutionState {
public:
  GenericSolutionState(ParameterList &parameter);
  GenericSolutionState(const GenericSolutionState &other);
  virtual ~GenericSolutionState();

  void setRVE() { m_hasRVEs = true; };
  virtual auto getCopy() -> std::shared_ptr<GenericSolutionState>;

  virtual auto getType() -> SolutionTypes {
    return SolutionTypes::GenericSolutionState;
  }

  void request_element_data_field(indexType elementId, indexType fieldId,
                                  indexType rows, indexType cols);
  auto get_element_data_field(indexType elementId, indexType fieldId)
      -> Types::MatrixXX<prec> &;
  void set_element_data_field(indexType elementId, indexType fieldId,
                 Types::MatrixXX<prec> &data);

  virtual void setValues(std::map<std::string, prec> &values){};

  virtual void setSparseMatrix(PointerCollection &pointers){};
  virtual void insertStiffnessResidual(
      Eigen::Matrix<prec, Eigen::Dynamic, Eigen::Dynamic> &stiffness,
      Eigen::Matrix<prec, Eigen::Dynamic, 1> &residual,
      std::vector<DegreeOfFreedom *> &Dofs);

  virtual void computeLoads(PointerCollection &pointers){};
  virtual void assembleSystem(PointerCollection &pointers){};
  virtual void factorize();
  virtual void solve(PointerCollection &pointers);
  virtual prec residual() { return this->residualVal; };
  virtual prec lgsResidual() { return 0; };
  virtual prec energyNorm() { return this->energyVal; };
  virtual void setEquationZero();
  virtual void updateSolution(PointerCollection &pointers){};
  virtual void dampedSolutionUpdate(PointerCollection &pointers){};
  virtual void computeEigenValues(PointerCollection &pointers, indexType number,
                                  indexType addNumber = 0, bool max = false,
                                  prec tol = 1e-10, prec shift = 1e-10){};

  virtual void setInitialValues(indexType numberOfEquations,
                                indexType numberOfActiveEquations);

  void setSolver(SolverTypes type);

  void setProps(std::shared_ptr<PropfunctionHandler> &props);
  std::shared_ptr<PropfunctionHandler> getProps();


  virtual auto getSolution(GenericNodes &node) -> Types::Vector3<prec> {
    return {0, 0, 0};
  };

  virtual prec getSolution(indexType globalId) { return 0; };
  virtual Eigen::Matrix<prec, Eigen::Dynamic, 1>
  getSolution(std::vector<DegreeOfFreedom *> &Dofs) {
    return Eigen::Matrix<prec, Eigen::Dynamic, 1>(0);
  };
  virtual Eigen::Matrix<prec, Eigen::Dynamic, 1>
  getIncrementalSolution(std::vector<DegreeOfFreedom *> &Dofs) {
    return Eigen::Matrix<prec, Eigen::Dynamic, 1>(0);
  };
  virtual Eigen::Matrix<prec, Eigen::Dynamic, 1>
  getNewtonSolution(std::vector<DegreeOfFreedom *> &Dofs) {
    return Eigen::Matrix<prec, Eigen::Dynamic, 1>(0);
  };

  virtual Eigen::Matrix<prec, Eigen::Dynamic, 1>
  getVelocity(std::vector<DegreeOfFreedom *> Dofs) {
    return Eigen::Matrix<prec, Eigen::Dynamic, 1>(0);
  };
  virtual Eigen::Matrix<prec, Eigen::Dynamic, 1>
  getAcceleration(std::vector<DegreeOfFreedom *> Dofs) {
    return Eigen::Matrix<prec, Eigen::Dynamic, 1>(0);
  };

  virtual void resetSolution() {
    this->eigenVectors.clear();
    this->eigenValues.clear();
  };

  virtual indexType numberOfEigenValues() {
    return static_cast<indexType>(this->eigenValues.size());
  };
  virtual prec getEigenVectorComp(indexType eqId, indexType evId);

  virtual void printSpMat(PointerCollection &pointers);

  virtual void computeConditionNumber(PointerCollection &pointers){};

  prec getEigenValue(indexType &number) {
    if (number < static_cast<indexType>(this->eigenValues.size())) {
      return this->eigenValues[number];
    } else {
      return (prec)0;
    }
  };

  virtual void nextSolutionStep();

  virtual void ctestout(PointerCollection &pointers){};

  // Homogenization part
  virtual void setStrains(PointerCollection& pointers, Types::VectorX<prec> &strains){};
  virtual void initHomogenization(PointerCollection &pointers,
                                  indexType homogenizationType,
                                  ParameterList &parameters){};
  virtual void computeAMatrix(PointerCollection &pointers){};
  virtual void homogenize(PointerCollection& pointers){};

  struct HomogenizedData {
    Types::VectorX<prec> sigma;
    Types::MatrixXX<prec> C;
  };
  virtual auto getHomogenizedData() -> HomogenizedData {
    return HomogenizedData();
  };

  // History Data
  void setupHistoryData(PointerCollection &pointers);
  auto getHistoryData(indexType elementId) -> HistoryDataIterator;

  virtual void toFile(PointerCollection &pointers, std::ofstream &out);
  virtual void fromFile(PointerCollection &pointers, std::ifstream &in);

  
  virtual void RVEDatatoFile(PointerCollection &pointers, std::ofstream &out);
  virtual void RVEDatafromFile(PointerCollection &pointers, std::ifstream &in);

  void insertRowInEqSystem(Types::SparseMatrix<prec, indexType> &matrix,
                           Types::SparseVector<prec, indexType> &row,
                           indexType rownumber);
  void insertColumnInEqSystem(Types::SparseMatrix<prec, indexType> &matrix,
                           Types::SparseVector<prec, indexType> &col,
                           indexType colnumber);

  auto getNewConstraint(ConstraintTypes constraintType)
  -> std::shared_ptr<BaseConstraint>;
  auto getConstraint(indexType constraintNumber)
      -> std::shared_ptr<BaseConstraint>;
  auto getNumberOfConstraints() -> indexType;

  auto getConstraintHandler() -> ConstraintHandler&;

  void setHasRVE() { m_hasRVEs = true; };
  void updateRVEHistory(PointerCollection &pointers);

protected:
  void assembleCsrMatrix(
      Eigen::SparseMatrix<prec, 0, indexType> &SpMat,
      Eigen::Matrix<prec, Eigen::Dynamic, Eigen::Dynamic> &stiffness,
      std::vector<DegreeOfFreedom *> &Dofs);
  void assembleCsrMatrix(Eigen::Matrix<prec, Eigen::Dynamic, 1> &toInsert,
                         Eigen::Matrix<prec, Eigen::Dynamic, 1> &toAdd,
                         std::vector<DegreeOfFreedom *> &Dofs);

  void compressVector(PointerCollection &pointers,
                      Eigen::Matrix<prec, Eigen::Dynamic, 1> &compVect,
                      Eigen::Matrix<prec, Eigen::Dynamic, 1> &uncompVect,
                      bool compr);

  std::shared_ptr<PropfunctionHandler> props;
  indexType numberOfEquations, NumberOfActiveEquations,
      NumberOfinActiveEquations;
  std::shared_ptr<GenericSolver> solver;
  bool symmetricSolver;
  bool upper;

  std::vector<std::vector<prec>> eigenVectors;
  std::vector<prec> eigenValues;

  ConstraintHandler m_constraints;

protected:
  prec residualVal, energyVal;

private:
  HistoryDataManager m_HistoryData;
  bool m_hasRVEs;
  ElementDataFields m_ElementData;
};

} /* namespace HierAMuS */
