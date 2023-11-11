// Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause

#include "solver/TransientSolutionNewmark.h"
#include "solver/GenericSolver.h"

#include "pointercollection/pointercollection.h"

//Equations
#include "EquationHandler.h"

#include "LoadList.h"
#include "PropfunctionHandler.h"

#include "finiteElements/ElementList.h"
#include "finiteElements/GenericFiniteElement.h"

#include "control/ParameterList.h"

namespace HierAMuS {

TransientSolutionNewmark::TransientSolutionNewmark(ParameterList &parameter)
    : TransientSolution(parameter) {

  if (parameter.hasParameter("beta")) {
    this->beta = parameter.getPrecVal("beta");
  } else {
    this->beta = prec(0.25);
  }
  if (parameter.hasParameter("gamma")) {
    this->gamma = parameter.getPrecVal("gamma");
  } else {
    this->gamma = prec(0.5);
  }

}

TransientSolutionNewmark::~TransientSolutionNewmark() {

}

void TransientSolutionNewmark::setValues(std::map<std::string, prec> &values) {
  typename std::map<std::string, prec>::iterator it;
  it = values.find("beta");
  it != values.end() ?
      this->beta = it->second :
      this->beta = (prec) 1 / (prec) 4;
  it = values.find("gamma");
  it != values.end() ?
      this->gamma = it->second :
      this->gamma = (prec) 1 / (prec) 2;
}

void TransientSolutionNewmark::setInitialValues(indexType numberOfEquations,
                                                indexType numberOfActiveEquations) {
  this->numberOfEquations = numberOfEquations;
  this->NumberOfActiveEquations = numberOfActiveEquations;
  this->Solution.resize(numberOfEquations);
  this->IncSolution.resize(numberOfEquations);
  this->dIncSolution.resize(numberOfEquations);
  this->vn.resize(numberOfEquations);
  this->vn1.resize(numberOfEquations);
  this->an.resize(numberOfEquations);
  this->an1.resize(numberOfEquations);

  this->SpMat.resize(this->NumberOfActiveEquations, this->NumberOfActiveEquations);
  this->Stiffness.resize(this->NumberOfActiveEquations, this->NumberOfActiveEquations);
  this->Mass.resize(this->NumberOfActiveEquations, this->NumberOfActiveEquations);
  this->Damping.resize(this->NumberOfActiveEquations, this->NumberOfActiveEquations);

  this->tripletList.reserve(numberOfActiveEquations);
  this->Rhs.resize(numberOfActiveEquations);
  this->eqSol.resize(numberOfActiveEquations);

  this->Solution.setZero();
  this->vn.setZero();
  this->vn1.setZero();
  this->an.setZero();
  this->an1.setZero();
}


void TransientSolutionNewmark::setSparseMatrix(PointerCollection& pointers) {
  GenericSolutionState::setSparseMatrix(pointers);
  this->SpMat.setFromTriplets(this->tripletList.begin(), this->tripletList.end());
  this->tripletList.clear();
  this->Mass = this->SpMat;
  this->Stiffness = this->SpMat;
  this->Damping = this->SpMat;
  this->solver->analyze(this->SpMat);
  std::cout << "Transient Algorthm: Newmark\n" << "Number of Equations: " << this->NumberOfActiveEquations << "\n" <<
            "Total number of degrees of Freedom: " << this->numberOfEquations << "\n" <<
            "Number of non-zero terms: " << this->SpMat.nonZeros() << std::endl;
}

void TransientSolutionNewmark::insertStiffnessResidual(Eigen::Matrix<prec, Eigen::Dynamic, Eigen::Dynamic> &stiffness,
                                                       Eigen::Matrix<prec, Eigen::Dynamic, 1> &residual,
                                                       std::vector<DegreeOfFreedom *> &Dofs) {

  indexType vsize = static_cast<indexType>(Dofs.size());

#pragma omp critical
  for (indexType i = 0; i < vsize; ++i) {
    if (Dofs[i]->getStatus() != dofStatus::inactive) {
      for (indexType j = 0; j < vsize; ++j) {
        if (Dofs[j]->getStatus() != dofStatus::inactive) {
          if (this->symmetricSolver) {
            if (this->upper) {
              if (Dofs[i]->getEqId() <= Dofs[j]->getEqId()) {
#pragma omp critical (assembleOne)
                this->Stiffness.coeffRef(Dofs[i]->getEqId(), Dofs[j]->getEqId()) += stiffness(i, j);
              }
            } else {
              if (Dofs[i]->getEqId() >= Dofs[j]->getEqId()) {
#pragma omp critical (assembleTwo)
                this->Stiffness.coeffRef(Dofs[i]->getEqId(), Dofs[j]->getEqId()) += stiffness(i, j);
              }
            }
          } else {
#pragma omp critical (assembleOne)
            this->Stiffness.coeffRef(Dofs[i]->getEqId(), Dofs[j]->getEqId()) += stiffness(i, j);
          }
        } else {
          this->Rhs(Dofs[i]->getEqId()) -= stiffness(i, j) * this->dIncSolution(Dofs[j]->getId());
        }
      }
#pragma omp critical (assembleThree)
      this->Rhs(Dofs[i]->getEqId()) -= residual(i);
    }
  }
  //Eigen::SparseMatrix<prec, 0, indexType> SpMat
}

void TransientSolutionNewmark::computeLoads(PointerCollection &pointers) {
  std::vector<prec> load, loadinc;
  std::vector<indexType> ids;

  auto eqHandler = pointers.getEquationHandler();

  pointers.getLoadList()->computeLoads(*this->getProps(), ids, load, loadinc);
  for (auto i = 0; i < ids.size(); ++i) {
    auto &dof = eqHandler->getDegreeOfFreedom(ids[i]);
    if (dof.getStatus() != dofStatus::inactive) {
      this->Rhs(dof.getEqId()) += load[i];
    } else {
      this->dIncSolution(dof.getId()) += loadinc[i];
    }
  }
}

void TransientSolutionNewmark::assembleSystem(PointerCollection& pointers) {
  indexType numberOfElements;
  auto elemList = pointers.getElementList();
  numberOfElements = elemList->getNumberOfElements();
  FiniteElement::GenericFiniteElement *elem;
  this->setEquationZero();

  Eigen::Matrix<prec, Eigen::Dynamic, Eigen::Dynamic> stiffness;
  Eigen::Matrix<prec, Eigen::Dynamic, 1> residual;
  std::vector<DegreeOfFreedom *> Dofs;

  this->computeLoads(pointers);

  //#pragma omp parallel for private(elem,stiffness,residual, Dofs) schedule(dynamic)
  for (indexType i = 0; i < numberOfElements; ++i) {
    Dofs.clear();
    elem = elemList->getElement(pointers, i);
    elem->GenericSetTangentResidual(pointers, stiffness, residual, Dofs);
    this->insertStiffnessResidual(stiffness, residual, Dofs);
    stiffness.setZero();
    residual.setZero();
    Dofs.clear();
    elem->GenericSetMass(pointers, stiffness, residual,Dofs);
    this->insertMassMatrix(stiffness, Dofs);
    //			//		elementLibrary(*elem, ControlOptions::BuildStiffnessResidual);
  }
  //
  prec dt = this->props->getDt();
  this->SpMat = (prec) 1 / this->beta / dt / dt * this->Mass + this->Stiffness;
  this->compressVector(pointers, this->compW, this->an1, true);
  //this->compW *= (prec) 1/beta/dt/dt;
  this->Rhs -= this->Mass * this->compW;
  this->residualVal = this->Rhs.norm();
}

void TransientSolutionNewmark::setEquationZero() {
  this->Stiffness.setZero();
  this->Damping.setZero();
  this->Mass.setZero();
  this->Rhs.setZero();
  this->dIncSolution.setZero();

}

Eigen::Matrix<prec, Eigen::Dynamic, 1> TransientSolutionNewmark::getSolution(std::vector<DegreeOfFreedom *> &Dofs) {
  Eigen::Matrix<prec, Eigen::Dynamic, 1> RetVec;

  indexType Size = static_cast<indexType>(Dofs.size());
  if (Size > 0) {
    RetVec.resize(Size);
    for (auto i = 0; i < Size; ++i) {
      RetVec(i) = this->Solution(Dofs[i]->getId());
    }
  }

  return RetVec;
}

Eigen::Matrix<prec, Eigen::Dynamic, 1> TransientSolutionNewmark::getVelocity(std::vector<DegreeOfFreedom *> Dofs) {
  Eigen::Matrix<prec, Eigen::Dynamic, 1> RetVec;

  indexType Size = static_cast<indexType>(Dofs.size());
  if (Size > 0) {
    RetVec.resize(Size);
    for (auto i = 0; i < Size; ++i) {
      RetVec(i) = this->vn(Dofs[i]->getId());
    }
  }

  return RetVec;
}

Eigen::Matrix<prec, Eigen::Dynamic, 1> TransientSolutionNewmark::getAcceleration(std::vector<DegreeOfFreedom *> Dofs) {
  Eigen::Matrix<prec, Eigen::Dynamic, 1> RetVec;

  indexType Size = static_cast<indexType>(Dofs.size());
  if (Size > 0) {
    RetVec.resize(Size);
    for (auto i = 0; i < Size; ++i) {
      RetVec(i) = this->an(Dofs[i]->getId());
    }
  }

  return RetVec;
}

void TransientSolutionNewmark::insertMassMatrix(Eigen::Matrix<prec, Eigen::Dynamic, Eigen::Dynamic> &stiffness,
                                                std::vector<DegreeOfFreedom *> &Dofs) {
  indexType vsize = static_cast<indexType>(Dofs.size());
#pragma omp critical
  for (indexType i = 0; i < vsize; ++i) {
    if (Dofs[i]->getStatus() != dofStatus::inactive) {
      for (indexType j = 0; j < vsize; ++j) {
        if (Dofs[j]->getStatus() != dofStatus::inactive) {
          if (this->symmetricSolver) {
            if (this->upper) {
              if (Dofs[i]->getEqId() <= Dofs[j]->getEqId()) {
#pragma omp critical (MassassembleOne)
                this->Mass.coeffRef(Dofs[i]->getEqId(), Dofs[j]->getEqId()) += stiffness(i, j);
              }
            } else {
              if (Dofs[i]->getEqId() >= Dofs[j]->getEqId()) {
#pragma omp critical (MassassembleTwo)
                this->Mass.coeffRef(Dofs[i]->getEqId(), Dofs[j]->getEqId()) += stiffness(i, j);
              }
            }
          } else {
#pragma omp critical (MassassembleOne)
            this->Mass.coeffRef(Dofs[i]->getEqId(), Dofs[j]->getEqId()) += stiffness(i, j);
          }
        }
      }
    }
  }
  //Eigen::SparseMatrix<prec, 0, indexType> SpMat
}

void TransientSolutionNewmark::nextSolutionStep() {
  this->props->incrTime();
  prec dt = this->props->getDt();
  this->an = this->an1;
  this->vn = this->vn1;
  prec one = static_cast<prec>(1);
  prec two = static_cast<prec>(2);
  this->vn1 = (one - this->gamma / this->beta) * this->vn + (one - this->gamma / two / this->beta) * this->an * dt;
  this->an1 = (one - one / two / this->beta) * this->an - one / this->beta / dt * this->vn;
}

void TransientSolutionNewmark::factorize() {
  this->solver->factorize(this->SpMat);
}

void TransientSolutionNewmark::solve(PointerCollection& pointers) {
  prec residual = this->Rhs.dot(this->Rhs);
  std::cout << "Residual norm: " << sqrt(residual) << std::endl;
  this->solver->solve(this->Rhs, this->eqSol);
  this->energyVal = sqrt(this->Rhs.dot(this->eqSol));
  this->props->update();
}

void HierAMuS::TransientSolutionNewmark::updateSolution(PointerCollection& pointers) {
  prec dt = this->props->getDt();
  this->compressVector(pointers, this->eqSol, this->dIncSolution, false);
  this->IncSolution += this->dIncSolution;
  this->Solution += this->dIncSolution;
  this->vn1 += this->gamma / this->beta / dt * this->dIncSolution;
  this->an1 += (prec) 1 / this->beta / dt / dt * this->dIncSolution;
}

prec TransientSolutionNewmark::getSolution(indexType globalId) {
  return this->Solution(globalId);
}

}


