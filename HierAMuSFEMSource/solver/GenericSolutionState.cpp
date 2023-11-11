// Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause

#include <pointercollection/pointercollection.h>

#include "finiteElements/GenericFiniteElement.h"


#include <control/ParameterList.h>
#include "PropfunctionHandler.h"

#include <solver/GenericSolutionState.h>

//Equations
#include "EquationHandler.h"

#include <solver/EigenPardisoLDLT.h>
#include <solver/EigenPardisoLLT.h>
#include <solver/EigenPardisoLU.h>
#include <solver/EigenSimplicialLDLT.h>
#include <solver/EigenSimplicialLLT.h>
#include <solver/EigenSparseLU.h>
#include <solver/EigenSparseQR.h>
#include <stdexcept>

#include "Constraints/GeneralLink.h"

#include "finiteElements/ElementList.h"

#include "control/BinaryWrite.h"

namespace HierAMuS {

GenericSolutionState::GenericSolutionState(ParameterList &parameter)
    : upper(false), symmetricSolver(false), solver(0),
      NumberOfActiveEquations(0), numberOfEquations(0), m_hasRVEs(false) {}

GenericSolutionState::GenericSolutionState(const GenericSolutionState &other)
{
  this->props = std::make_shared<PropfunctionHandler>((*other.props));
  this->numberOfEquations = other.numberOfEquations;
  this->NumberOfActiveEquations = other.NumberOfActiveEquations;
  this->NumberOfinActiveEquations = other.NumberOfinActiveEquations;
  this->setSolver(other.solver->getType());
  this->symmetricSolver = other.symmetricSolver;
  this->upper = other.upper;
  this->eigenVectors = other.eigenVectors;
  this->eigenValues = other.eigenValues;
  this->residualVal = other.residualVal;
  this->energyVal = other.energyVal;
  this->m_constraints = ConstraintHandler(other.m_constraints);
  this->m_HistoryData = other.m_HistoryData;
  this->m_hasRVEs = other.m_hasRVEs;
}

GenericSolutionState::~GenericSolutionState() {}

auto GenericSolutionState::getCopy() -> std::shared_ptr<GenericSolutionState> {
  std::shared_ptr<GenericSolutionState> tor =
      std::make_shared<GenericSolutionState>(*this);
  return tor;
}

void GenericSolutionState::request_element_data_field(indexType elementId,
                                                      indexType fieldId,
                                                      indexType rows,
                                                      indexType cols) {
  m_ElementData.request_field(elementId, fieldId, rows, cols);
}

auto GenericSolutionState::get_element_data_field(indexType elementId, indexType fieldId)
    -> Types::MatrixXX<prec> & {
  return m_ElementData.get_field(elementId, fieldId);
}

void GenericSolutionState::set_element_data_field(indexType elementId, indexType fieldId,
                                     Types::MatrixXX<prec> &data) {
  m_ElementData.set_field(elementId, fieldId, data);
}

void GenericSolutionState::setInitialValues(indexType numberOfEquations,
                                            indexType numberOfActiveEquations) {
}

void GenericSolutionState::insertStiffnessResidual(
    Eigen::Matrix<prec, Eigen::Dynamic, Eigen::Dynamic> &stiffness,
    Eigen::Matrix<prec, Eigen::Dynamic, 1> &residual,
    std::vector<DegreeOfFreedom *> &Dofs) {}

void GenericSolutionState::setEquationZero() {}

void GenericSolutionState::factorize() {}

void GenericSolutionState::solve(PointerCollection &pointers) {}

void GenericSolutionState::setSolver(SolverTypes type) {

  switch (type) {
  case SolverTypes::TypeEigenPardisoLDLT:
    this->solver = std::make_shared<EigenPardisoLDLT>();
    this->symmetricSolver = true;
    this->upper = true;
    break;
  case SolverTypes::TypeEigenPardisoLLT:
    this->solver = std::make_shared<EigenPardisoLLT>();
    this->symmetricSolver = true;
    this->upper = true;
    break;
  case SolverTypes::TypeEigenPardisoLU:
    this->solver = std::make_shared<EigenPardisoLU>();
    this->symmetricSolver = false;
    break;
  case SolverTypes::TypeEigenSimplicialLDLT:
    this->solver = std::make_shared<EigenSimplicialLDLT>();
    this->symmetricSolver = true;
    this->upper = false;
    break;
  case SolverTypes::TypeEigenSimplicialLLT:
    this->solver = std::make_shared<EigenSimplicialLLT>();
    this->symmetricSolver = true;
    this->upper = false;
    break;
  case SolverTypes::TypeEigenSparseLU:
    this->solver = std::make_shared<EigenSparseLU>();
    this->symmetricSolver = false;
    break;
  case SolverTypes::TypeEigenSparseQR:
    this->solver = std::make_shared<EigenSparseQR>();
    break;
  default:
    throw std::runtime_error("Solver type not implemented");
  }
}

void GenericSolutionState::setProps(
    std::shared_ptr<PropfunctionHandler> &props) {
  if (!this->props) {
    this->props = props;
  }
}

prec GenericSolutionState::getEigenVectorComp(indexType eqId, indexType evId) {
  if (eqId >= 0) {
    if (static_cast<indexType>(this->eigenVectors.size()) > evId) {
      return this->eigenVectors[evId][eqId];
    }
  }
  return prec(0);
}

void GenericSolutionState::printSpMat(PointerCollection &pointers) {

  // for (indexType k = 0; k < this->SpMat.outerSize(); ++k)
  //	for (Eigen::SparseMatrix<prec>::InnerIterator it(this->SpMat, k); it;
  //++it)
  //	{
  //		std::cout << it.row() << " " << it.col() << " " << it.value() <<
  // std::endl; 		it.value(); 		it.row();   // row index
  // it.col();
  // // col index (here it is equal to k) 		it.index(); // inner
  // index, here it is equal to it.row()
  //	}
}

void GenericSolutionState::setupHistoryData(PointerCollection &pointers) {
  auto elemList = pointers.getElementList();
  indexType numElem = elemList->getNumberOfElements();
  m_HistoryData.setNumberOfElements(numElem);
  for (indexType i = 0; i < numElem; i++) {
    auto elem = elemList->getElement(pointers, i);
    auto setupHist = elem->getHistoryDataSetUp(pointers);
    m_HistoryData.addHistoryDataElement(setupHist);
  }
  m_HistoryData.initHistoryData();
}

auto GenericSolutionState::getHistoryData(indexType elementId)
    -> HistoryDataIterator {
  return m_HistoryData.getHistoryDataIterator(elementId);
}

void GenericSolutionState::toFile(PointerCollection &pointers,
                                  std::ofstream &out) {
  // General Data
  writeScalar(out, numberOfEquations);
  writeScalar(out, NumberOfActiveEquations);
  writeScalar(out, NumberOfinActiveEquations);
  writeScalar(out, symmetricSolver);
  writeScalar(out, upper);
  writeScalar(out, m_hasRVEs);
  prec ctime = this->props->getTime();
  writeScalar(out, ctime);
  auto solverType = this->solver->getType();
  writeScalar(out, solverType);
  this->m_HistoryData.toFile(out);
  m_constraints.toFile(out);
  
}

void GenericSolutionState::fromFile(PointerCollection &pointers,
                                    std::ifstream &in) {
  // General Data
  readScalar(in, numberOfEquations);
  readScalar(in, NumberOfActiveEquations);
  readScalar(in, NumberOfinActiveEquations);
  readScalar(in, symmetricSolver);
  readScalar(in, upper);
  readScalar(in, m_hasRVEs);
  prec ctime;
  readScalar(in, ctime);
  SolverTypes solverType;
  readScalar(in, solverType);
  this->setSolver(solverType);

  this->props->setTime(ctime);
  this->m_HistoryData.fromFile(in);
  m_constraints.fromFile(in);

 
}

void GenericSolutionState::RVEDatatoFile(PointerCollection &pointers,
                                         std::ofstream &out) {
  prec ctime = this->props->getTime();
  writeScalar(out, ctime);
  auto solverType = this->solver->getType();
  writeScalar(out, solverType);
  this->m_HistoryData.toFile(out);
  //m_constraints.toFile(out);
}

void GenericSolutionState::RVEDatafromFile(PointerCollection &pointers,
                                           std::ifstream &in) {
  prec ctime;
  readScalar(in, ctime);
  SolverTypes solverType;
  readScalar(in, solverType);
  this->setSolver(solverType);

  this->props->setTime(ctime);
  this->m_HistoryData.fromFile(in);
  //m_constraints.fromFile(in);
}

void GenericSolutionState::insertRowInEqSystem(
    Types::SparseMatrix<prec, indexType> &matrix,
    Types::SparseVector<prec, indexType> &row, indexType rownumber) {
  if (this->symmetricSolver) {
    if (this->upper) {
      for (Types::SparseVector<prec, indexType>::InnerIterator it(row); it;
           ++it) {
        if (it.index() >= rownumber)
          matrix.coeffRef(rownumber, it.index()) += it.value();
      }
    } else {
      for (Types::SparseVector<prec, indexType>::InnerIterator it(row); it;
           ++it) {
        if (it.index() <= rownumber)
          matrix.coeffRef(rownumber, it.index()) += it.value();
      }
    }
  } else {
    for (Types::SparseVector<prec, indexType>::InnerIterator it(row); it;
         ++it) {
      matrix.coeffRef(rownumber, it.index()) += it.value();
    }
  }
}

void GenericSolutionState::insertColumnInEqSystem(
    Types::SparseMatrix<prec, indexType> &matrix,
    Types::SparseVector<prec, indexType> &col, indexType colnumber) {
  if (this->symmetricSolver) {
    if (this->upper) {
      for (Types::SparseVector<prec, indexType>::InnerIterator it(col); it;
           ++it) {
        if (it.index() <= colnumber)
          matrix.coeffRef(it.index(), colnumber) += it.value();
      }
    } else {

      for (Types::SparseVector<prec, indexType>::InnerIterator it(col); it;
           ++it) {
        if (it.index() >= colnumber)
          matrix.coeffRef(it.index(), colnumber) += it.value();
      }
    }
  } else {
    for (Types::SparseVector<prec, indexType>::InnerIterator it(col); it;
         ++it) {
      matrix.coeffRef(it.index(), colnumber) += it.value();
    }
  }
}

void GenericSolutionState::assembleCsrMatrix(
    Eigen::SparseMatrix<prec, 0, indexType> &SpMat,
    Eigen::Matrix<prec, Eigen::Dynamic, Eigen::Dynamic> &stiffness,
    std::vector<DegreeOfFreedom *> &Dofs) {
  DegreeOfFreedom *temp, *temp2;
  indexType vsize = static_cast<indexType>(Dofs.size());
  for (auto i = 0; i < vsize; ++i) {
    temp = Dofs[i];
    if (temp->getStatus() != dofStatus::inactive) {
      for (auto j = 0; j < vsize; ++j) {
        temp2 = Dofs[j];
        if (temp2->getStatus() != dofStatus::inactive) {
          bool add = false;
          indexType li = temp->getEqId(), lj = temp2->getEqId();

          this->symmetricSolver ? this->upper
                                      ? li <= lj ? add = true : add = false
                                  : li >= lj ? add = true
                                             : add = false
                                : add = true;

          if (add)
#pragma omp atomic
            SpMat.coeffRef(li, lj) += stiffness(i, j);
        }
      }
    }
  }
}

void HierAMuS::GenericSolutionState::assembleCsrMatrix(
    Eigen::Matrix<prec, Eigen::Dynamic, 1> &toInsert,
    Eigen::Matrix<prec, Eigen::Dynamic, 1> &residual,
    std::vector<DegreeOfFreedom *> &Dofs) {}


void GenericSolutionState::compressVector(
    PointerCollection &pointers,
    Eigen::Matrix<prec, Eigen::Dynamic, 1> &compVect,
    Eigen::Matrix<prec, Eigen::Dynamic, 1> &uncompVect, bool compr) {

  auto eqHandler = pointers.getEquationHandler();
  indexType numNodes = eqHandler->getNumberOfNodes();
  if (compr) {
    compVect.resize(this->NumberOfActiveEquations);
    for (auto i = 0; i < numNodes; ++i) {
      auto &node = eqHandler->getNode(i);
      for (auto j = 0; j < 3; ++j) {
        auto &dof = node.getDegreeOfFreedom(j);
        if (dof.getStatus() != dofStatus::inactive) {
          compVect(dof.getEqId()) = uncompVect(dof.getId());
        }
      }
    }
  } else {
    uncompVect.resize(this->numberOfEquations);
    for (auto i = 0; i < numNodes; ++i) {
      auto &node = eqHandler->getNode(i);
      for (auto j = 0; j < 3; ++j) {
        auto &dof = node.getDegreeOfFreedom(j);
        if (dof.getStatus() != dofStatus::inactive) {
          uncompVect(dof.getId()) = compVect(dof.getEqId());
        }
      }
    }
  }
}

std::shared_ptr<PropfunctionHandler> GenericSolutionState::getProps() {
  return this->props;
}

auto GenericSolutionState::getNewConstraint(ConstraintTypes constraintType)
-> std::shared_ptr<BaseConstraint> {
  
  return m_constraints.getNewConstraint(constraintType);
}

auto GenericSolutionState::getConstraint(indexType constraintNumber)
    -> std::shared_ptr<BaseConstraint> {
  return m_constraints.getConstraint(constraintNumber);
}

auto GenericSolutionState::getNumberOfConstraints() -> indexType {
  return m_constraints.getNumberOfConstraints();
}

auto GenericSolutionState::getConstraintHandler()
-> ConstraintHandler&
{
  return m_constraints;
}

void GenericSolutionState::updateRVEHistory(PointerCollection &pointers)
{

  if (m_hasRVEs) {
    indexType numberOfElements;
    auto elemList = pointers.getElementList();
    numberOfElements = elemList->getNumberOfElements();
    FiniteElement::GenericFiniteElement *elem;
    

    
#pragma omp parallel for               \
    schedule(dynamic)
    for (indexType i = 0; i < numberOfElements; ++i) {
      elem = elemList->getElement(pointers, i);
      elem->updateRVEHistory(pointers);
    }
  }
}

void GenericSolutionState::nextSolutionStep() {
  this->props->incrTime();
  m_HistoryData.update();

};

void ElementDataFields::request_field(indexType elementId, indexType fieldId,
                                      indexType rows, indexType cols) {
  indexType numField = m_data[elementId].size();
  if (fieldId + 1 > numField)
    m_data[elementId].resize(fieldId + 1);
  m_data[elementId][fieldId].resize(rows, cols);
  m_data[elementId][fieldId].setZero();
}

auto ElementDataFields::get_field(indexType elementId, indexType fieldId)
    -> Types::MatrixXX<prec> & {
  indexType numfields = m_data[elementId].size();
  if (fieldId + 1 > numfields)
    throw std::runtime_error("Error, requested field does not exist");
  return m_data[elementId][fieldId];
}

void ElementDataFields::set_field(indexType elementId, indexType fieldId,
                                  Types::MatrixXX<prec> &data) {
  indexType numFields = m_data[elementId].size();
  if (fieldId + 1 > numFields)
    m_data[elementId].resize(fieldId + 1);
  m_data[elementId][fieldId] = data;
}

} /* namespace HierAMuS */
