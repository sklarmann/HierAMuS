// Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause


#include "control/ParameterList.h"

#include <types/MatrixTypes.h>

#include <control/HandlingStructs.h>
#include <control/OutputHandler.h>
#include <pointercollection/pointercollection.h>

#include <solver/GenericSolutionState.h>
#include <solver/StaticSolutionState.h>

#include <control/HandlingStructs.h>
#include <pointercollection/pointercollection.h>

#include <finiteElements/ElementList.h>
#include "finiteElements/GenericFiniteElement.h"

#include <fstream>

//Equations
#include "EquationHandler.h"

#include "LoadList.h"
#include "PropfunctionHandler.h"

#include <solver/GenericSolver.h>

#include <sstream>
#include <vector>

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/SparseCore>

#include <sstream>

#include <cmath>

#include <Timer.h>

#include "spdlog/stopwatch.h"

#include "control/BinaryWrite.h"

#ifdef USE_SPECTRA
#include <Spectra/GenEigsRealShiftSolver.h>
#include <Spectra/GenEigsSolver.h>
#include <Spectra/MatOp/SparseGenMatProd.h>
#include <Spectra/MatOp/SparseGenRealShiftSolve.h>
#include <Spectra/MatOp/SparseSymMatProd.h>
#include <Spectra/MatOp/SparseSymShiftSolve.h>
#include <Spectra/MatOp/SparseGenMatProd.h>
#include <Spectra/MatOp/SparseGenRealShiftSolve.h>
#include <Spectra/SymEigsShiftSolver.h>
#include <Spectra/SymEigsSolver.h>
#endif // USE_SPECTRA

#include <iostream>

namespace HierAMuS {

StaticSolutionState::StaticSolutionState(ParameterList &parameter)
    : GenericSolutionState(parameter) {}

StaticSolutionState::StaticSolutionState(const StaticSolutionState &other)
  : GenericSolutionState(other)
{
  this->Solution = other.Solution;
  this->IncSolution = other.IncSolution;
  this->dIncSolution = other.dIncSolution;
  this->NewtonSolution = other.NewtonSolution;

  this->Rhs = other.Rhs;
  this->RhsB = other.RhsB;
  this->eqSol = other.eqSol;
  this->incSol = other.incSol;

  this->Kaa = other.Kaa;
  this->Kab = other.Kab;
  this->Kbb = other.Kbb;
  this->Kba = other.Kba;

  this->solver->analyze(this->Kaa);
}

StaticSolutionState::~StaticSolutionState() { this->Kaa.setZero(); }

auto StaticSolutionState::getCopy() -> std::shared_ptr<GenericSolutionState> {

  std::shared_ptr<GenericSolutionState> tor =
      std::make_shared<StaticSolutionState>(*this);
  return tor;
}

void StaticSolutionState::setInitialValues(indexType numberOfEquations,
                                           indexType numberOfActiveEquations) {
  this->numberOfEquations = numberOfEquations;
  this->NumberOfActiveEquations = numberOfActiveEquations;
  this->NumberOfinActiveEquations = numberOfEquations - numberOfActiveEquations;

  this->Solution.resize(numberOfEquations);
  this->IncSolution.resize(numberOfEquations);
  this->dIncSolution.resize(numberOfEquations);
  this->NewtonSolution.resize(numberOfEquations);
  this->Kaa.resize(this->NumberOfActiveEquations,
                   this->NumberOfActiveEquations);
  this->Kab.resize(this->NumberOfActiveEquations,
                   this->NumberOfinActiveEquations);
  this->Kba.resize(this->NumberOfinActiveEquations,
                   this->NumberOfActiveEquations);
  this->Kbb.resize(this->NumberOfinActiveEquations,
                   this->NumberOfinActiveEquations);

  this->Rhs.resize(numberOfActiveEquations);
  this->RhsB.resize(NumberOfinActiveEquations);
  this->eqSol.resize(numberOfActiveEquations);
  this->incSol.resize(NumberOfinActiveEquations);

  this->Solution.setZero();
  this->NewtonSolution.setZero();
}

void StaticSolutionState::setSparseMatrix(PointerCollection &pointers) {

  Types::VectorX<indexType> entriesPerColActive, entriesPerColinActive,
      entriesPerColLowerLeft;
  auto elemList = pointers.getElementList();
  auto numElem = elemList->getNumberOfElements();
  auto equ = pointers.getEquationHandler();
  equ->updateEquations();
  this->setInitialValues(equ->getNumberOfTotalEquations(),
                         equ->getNumberOfActiveEquations());
 
  m_constraints.initialize(pointers);

  entriesPerColActive.resize(this->NumberOfActiveEquations);
  entriesPerColActive.setZero();

  entriesPerColLowerLeft.resize(this->NumberOfActiveEquations);
  entriesPerColLowerLeft.setZero();

  entriesPerColinActive.resize(this->NumberOfinActiveEquations);
  entriesPerColinActive.setZero();

  for (auto i = 0; i < numElem; ++i) {
    auto elem = elemList->getElement(pointers, i);
    auto Dofs = elem->getDofs(pointers);
    indexType activeCols = 0;
    indexType inactiveCols = 0;
    for (auto j : Dofs) {
      if (j->getStatus() == dofStatus::active) {
        activeCols++;
      } else {
        inactiveCols++;
      }
    }

    for (auto j : Dofs) {
      if (j->getStatus() == dofStatus::active) {
        entriesPerColActive(j->getEqId()) += activeCols;
        entriesPerColLowerLeft(j->getEqId()) += inactiveCols;
      } else {
        entriesPerColinActive(j->getEqId()) += inactiveCols;
      }
    }
  }
  

  this->Kaa.reserve(entriesPerColActive);
  this->Kba.reserve(entriesPerColLowerLeft);
  this->Kbb.reserve(entriesPerColinActive);

  for (auto i = 0; i < numElem; ++i) {
    auto Dofs = elemList->getElement(pointers, i)->getDofs(pointers);
    for (auto j : Dofs) {
      if (j->getStatus() == dofStatus::active) { //(Kaa) j & k are active
        for (auto k : Dofs) {
          if (k->getStatus() == dofStatus::active) {
            if (this->symmetricSolver) {
              if (this->upper) {
                if (j->getEqId() <= k->getEqId())
                  this->Kaa.coeffRef(j->getEqId(), k->getEqId()) = 1.0;
              } else {
                if (j->getEqId() >= k->getEqId())
                  this->Kaa.coeffRef(j->getEqId(), k->getEqId()) = 1.0;
              }
            } else {
              this->Kaa.coeffRef(j->getEqId(), k->getEqId()) = 1.0;
            }
          }
        }
      } else {
        //(Kba) j not active & k is active
        for (auto k : Dofs) {
          if (k->getStatus() == dofStatus::active) {
            this->Kba.coeffRef(j->getEqId(), k->getEqId()) = 1.0;
          } else {
            if (this->symmetricSolver) {
              if (this->upper) {
                if (j->getEqId() <= k->getEqId())
                  this->Kbb.coeffRef(j->getEqId(), k->getEqId()) = 1.0;
              } else {
                if (j->getEqId() >= k->getEqId())
                  this->Kbb.coeffRef(j->getEqId(), k->getEqId()) = 1.0;
              }
            } else {
              this->Kbb.coeffRef(j->getEqId(), k->getEqId()) = 1.0;
            }
          }
        }
      }
    }
  }

  this->Kab = this->Kba.transpose();
                                                             
  m_constraints.modifyEquationSystem(pointers, Kaa, Kab, Kba, Kbb, Rhs, RhsB,
                                     symmetricSolver, upper);

  this->Kaa.makeCompressed();
  this->solver->analyze(this->Kaa);
  this->Kab.makeCompressed();
  this->Kba.makeCompressed();
  this->Kbb.makeCompressed();
}

void StaticSolutionState::insertStiffnessResidual(
    Eigen::Matrix<prec, Eigen::Dynamic, Eigen::Dynamic> &stiffness,
    Eigen::Matrix<prec, Eigen::Dynamic, 1> &residual,
    std::vector<DegreeOfFreedom *> &Dofs) {

  if (this->symmetricSolver)
  {
    if (this->upper)
    {
      insertStiffnessResidualSymUpper(stiffness, residual, Dofs);
    } else
    {
      insertStiffnessResidualSymLower(stiffness, residual, Dofs);
    }
  } else
  {
    insertStiffnessResidualFull(stiffness,residual,Dofs);
  }
}

void StaticSolutionState::insertStiffnessResidualSymUpper(
    Eigen::Matrix<prec, Eigen::Dynamic, Eigen::Dynamic> &stiffness,
    Eigen::Matrix<prec, Eigen::Dynamic, 1> &residual,
    std::vector<DegreeOfFreedom *> &Dofs) {
  auto vsize = Dofs.size();

  for (auto i = 0; i < vsize; ++i) {
    auto eqI = static_cast<Eigen::Index>(Dofs[i]->getEqId());
    if (Dofs[i]->getStatus() == dofStatus::active) {
      for (auto j = 0; j < vsize; ++j) {
        auto eqJ = static_cast<Eigen::Index>(Dofs[j]->getEqId());
        if (Dofs[j]->getStatus() == dofStatus::active) {
          if (eqI <= eqJ) {
#pragma omp atomic
            this->Kaa.coeffRef(eqI, eqJ) += stiffness(i, j);
          }
        } // end active j
        else {
#pragma omp atomic
          this->Kab.coeffRef(eqI, eqJ) += stiffness(i, j);
//#pragma omp atomic
//          this->Rhs(eqI) -= stiffness(i,j)*this->incSol(Dofs[j]->getEqId());
        }
      }   // end j loop
#pragma omp atomic
      this->Rhs(eqI) -= residual(i);
    } // end active i
    else {
      for (auto j = 0; j < vsize; ++j) {
        auto eqJ = static_cast<Eigen::Index>(Dofs[j]->getEqId());
        if (Dofs[j]->getStatus() == dofStatus::active) {
#pragma omp atomic
          Kba.coeffRef(eqI, eqJ) += stiffness(i, j);
        } else {
          if (eqI <= eqJ) {
#pragma omp atomic
            Kbb.coeffRef(eqI, eqJ) += stiffness(i, j);
          }
        }
      }
#pragma omp atomic
      this->RhsB(eqI) -= residual(i);
    }
  } // end i loop
}
void StaticSolutionState::insertStiffnessResidualSymLower(
    Eigen::Matrix<prec, Eigen::Dynamic, Eigen::Dynamic> &stiffness,
    Eigen::Matrix<prec, Eigen::Dynamic, 1> &residual,
    std::vector<DegreeOfFreedom *> &Dofs) {
  indexType vsize = static_cast<indexType>(Dofs.size());

  for (indexType i = 0; i < vsize; ++i) {
    indexType eqI = Dofs[i]->getEqId();
    if (Dofs[i]->getStatus() == dofStatus::active) {
      for (indexType j = 0; j < vsize; ++j) {
        indexType eqJ = Dofs[j]->getEqId();
        if (Dofs[j]->getStatus() == dofStatus::active) {
          if (eqI >= eqJ) {
#pragma omp atomic
            this->Kaa.coeffRef(eqI, eqJ) += stiffness(i, j);
          }
        } // end active j
        else {
#pragma omp atomic
          this->Kab.coeffRef(eqI, eqJ) += stiffness(i, j);
//#pragma omp atomic
//          this->Rhs(eqI) -= stiffness(i,j)*this->incSol(Dofs[j]->getEqId());
        }
      }   // end j loop
#pragma omp atomic
      this->Rhs(eqI) -= residual(i);
    } // end active i
    else {
      for (indexType j = 0; j < vsize; ++j) {
        indexType eqJ = Dofs[j]->getEqId();
        if (Dofs[j]->getStatus() == dofStatus::active) {
#pragma omp atomic
          Kba.coeffRef(eqI, eqJ) += stiffness(i, j);
        } else {
          if (eqI >= eqJ) {
#pragma omp atomic
            Kbb.coeffRef(eqI, eqJ) += stiffness(i, j);
          }
        }
      }
#pragma omp atomic
      this->RhsB(eqI) -= residual(i);
    }
  } // end i loop
}
void StaticSolutionState::insertStiffnessResidualFull(
    Eigen::Matrix<prec, Eigen::Dynamic, Eigen::Dynamic> &stiffness,
    Eigen::Matrix<prec, Eigen::Dynamic, 1> &residual,
    std::vector<DegreeOfFreedom *> &Dofs) {
  indexType vsize = static_cast<indexType>(Dofs.size());

  for (indexType i = 0; i < vsize; ++i) {
    indexType eqI = Dofs[i]->getEqId();
    if (Dofs[i]->getStatus() == dofStatus::active) {
      for (indexType j = 0; j < vsize; ++j) {
        indexType eqJ = Dofs[j]->getEqId();
        if (Dofs[j]->getStatus() == dofStatus::active) {
#pragma omp atomic
          this->Kaa.coeffRef(eqI, eqJ) += stiffness(i, j);

        } // end active j
        else 
        {
#pragma omp atomic
          this->Kab.coeffRef(eqI, eqJ) += stiffness(i, j);
        }
      }   // end j loop
#pragma omp atomic
      this->Rhs(eqI) -= residual(i);
    } // end active i
    else {
      for (indexType j = 0; j < vsize; ++j) {
        indexType eqJ = Dofs[j]->getEqId();
        if (Dofs[j]->getStatus() == dofStatus::active) {
#pragma omp atomic
          Kba.coeffRef(eqI, eqJ) += stiffness(i, j);
        } else {
#pragma omp atomic
          Kbb.coeffRef(eqI, eqJ) += stiffness(i, j);
        }
      }
#pragma omp atomic
      this->RhsB(eqI) -= residual(i);
    }
  } // end i loop
}


void StaticSolutionState::setEquationZero() {
  auto SpMatZero = [](Types::SparseMatrix<prec, indexType> &matrix) {
    typedef
        typename Eigen::SparseMatrix<prec, 0, indexType>::InnerIterator eigenIt;
    for (auto k = 0; k < matrix.outerSize(); ++k) {
      for (eigenIt it(matrix, k); it; ++it) {
        it.valueRef() = 0.0;
      }
    }
  };
  SpMatZero(Kaa);
  SpMatZero(Kab);
  SpMatZero(Kba);
  SpMatZero(Kbb);

  RhsB.setZero();
  this->Rhs.setZero();
  this->dIncSolution.setZero();
  this->incSol.setZero();
}

void StaticSolutionState::factorize() { this->solver->factorize(this->Kaa); }

void StaticSolutionState::solve(PointerCollection &pointers) {

  this->solver->solve(this->Rhs, this->eqSol);
  prec dotpp = this->Rhs.dot(this->eqSol);
  dotpp = abs(dotpp);
  this->energyVal = sqrt(dotpp);
  pointers.getPropLoads()->update();
}

auto StaticSolutionState::lgsResidual() -> prec {
  prec temp;
  if (this->symmetricSolver) {
    if (this->upper) {
      temp = (this->Kaa.template selfadjointView<Eigen::Upper>() * this->eqSol -
              this->Rhs)
                 .norm();
    } else {
      temp = (this->Kaa.template selfadjointView<Eigen::Lower>() * this->eqSol -
              this->Rhs)
                 .norm();
      //((this->Kaa).selfadjointView<Eigen::Lower>());
    }
  } else {
    temp = (this->Kaa * this->eqSol - this->Rhs).norm();
  }
  return temp;
}

void StaticSolutionState::computeLoads(PointerCollection &pointers) {
  std::vector<prec> load, loadinc;
  std::vector<indexType> ids;

  auto eqHandler = pointers.getEquationHandler();

  pointers.getLoadList()->computeLoads(*this->getProps(), ids, load, loadinc);
  for (auto i = 0; i < ids.size(); ++i) {
    auto &dof = eqHandler->getDegreeOfFreedom(ids[i]);
    if (dof.getStatus() == dofStatus::active) {
      this->Rhs(dof.getEqId()) += load[i];
    } else {
      this->RhsB(dof.getEqId()) += load[i];
    }
  }
  load.clear();
  loadinc.clear();
  ids.clear();
  pointers.getPrescribedDisplacements()->computeLoads(*this->getProps(), ids, load,
                                                      loadinc);
  for (auto i = 0; i < ids.size(); ++i) {
    auto &dof = eqHandler->getDegreeOfFreedom(ids[i]);
    if (dof.getStatus() == dofStatus::active) {
      pointers.getSPDLogger().warn("In StaticSolutionState method computeLoads, "
                                   "something wrong with prescribed "
                                   "displacements. Dof status of dof with id: "
                                   "{} is active instead of inactive.",
                                   dof.getId());
    } else {
      //this->dIncSolution(dof.getId()) += loadinc[i];
      this->incSol(dof.getEqId()) += load[i]-this->getSolution(dof.getId());
    }
  }
}

void StaticSolutionState::assembleSystem(PointerCollection &pointers) {

  indexType numberOfElements;
  auto elemList = pointers.getElementList();
  numberOfElements = elemList->getNumberOfElements();
  FiniteElement::GenericFiniteElement *elem;
  this->setEquationZero();

  Eigen::Matrix<prec, Eigen::Dynamic, Eigen::Dynamic> stiffness;
  Eigen::Matrix<prec, Eigen::Dynamic, 1> residual;
  std::vector<DegreeOfFreedom *> Dofs;

  this->computeLoads(pointers);

#pragma omp parallel for private(elem, stiffness, residual, Dofs)              \
    schedule(dynamic)
  for (indexType i = 0; i < numberOfElements; ++i) {
    Dofs.clear();
    elem = elemList->getElement(pointers, i);
    elem->GenericSetTangentResidual(pointers, stiffness, residual, Dofs);
    this->insertStiffnessResidual(stiffness, residual, Dofs);
  }

  Timer<sec> ttimer;
  spdlog::stopwatch ws;    
  
  ttimer.start();
  m_constraints.modifyEquationSystem(pointers, Kaa, Kab, Kba, Kbb, Rhs, RhsB,
                                     symmetricSolver, upper);
  ttimer.stop();

  this->Rhs -= this->Kab * this->incSol;

  pointers.getSPDLogger().debug("    Modification of eq system took:    {} sec",ws);
  this->residualVal = this->Rhs.norm();
}

void StaticSolutionState::nextSolutionStep() {
  GenericSolutionState::nextSolutionStep();
  this->IncSolution.setZero();
  this->NewtonSolution.setZero();
  this->dIncSolution.setZero();
}

void StaticSolutionState::updateSolution(PointerCollection &pointers) {
  auto eqHandler = pointers.getEquationHandler();
  indexType numberOfNodes = eqHandler->getNumberOfNodes();

  for (auto i = 0; i < numberOfNodes; ++i) {
    auto &node = eqHandler->getNode(i);
    for (auto j = 0; j < 3; ++j) {
      auto &dof = node.getDegreeOfFreedom(j);
      indexType id = dof.getId();
      indexType eqId = dof.getEqId();
      if (dof.getStatus() == dofStatus::active) {
        this->NewtonSolution(id) = this->eqSol(eqId);
      } else {
        this->NewtonSolution(id) = this->incSol(eqId);
      }
    }
  }
  
  this->Solution += this->NewtonSolution;
  this->NewtonSolution = m_constraints.getSlaveNewtonSolution(pointers);
  
  this->Solution += this->NewtonSolution;
  this->IncSolution += this->NewtonSolution;
}

void StaticSolutionState::dampedSolutionUpdate(PointerCollection &pointers) {
  prec cres = this->residual();
  Eigen::Matrix<prec, Eigen::Dynamic, 1> BSolution, BIncSolution, BdIncSolution,
      BNewtonSolution;
  BSolution = this->Solution;
  BIncSolution = this->IncSolution;
  BdIncSolution = this->dIncSolution;
  BNewtonSolution = this->NewtonSolution;
  this->updateSolution(pointers);
  this->assembleSystem(pointers);
  prec nres = this->residual();
  // indexType cc = 0;
  while (nres > cres) {
    // cc += 1;
    this->Solution = BSolution;
    this->IncSolution = BIncSolution;
    this->dIncSolution = BdIncSolution;
    this->NewtonSolution = BNewtonSolution;
    this->eqSol /= 2;
    this->updateSolution(pointers);
    this->assembleSystem(pointers);
    nres = this->residual();
  }
}

void StaticSolutionState::computeEigenValues(PointerCollection &pointers,
                                             indexType number,
                                             indexType addNumber /* = 0 */,
                                             bool max /* = false */,
                                             prec tol /* = 1e-10 */
                                             ,
                                             prec shift /* = 1e-10 */) {

  if (this->symmetricSolver)
  {
    if (this->upper) {
      this->computeEigenValuesSymUpper(pointers, number, addNumber, max, tol,
                                       shift);
      
    }else
    {
      this->computeEigenValuesSymLower(pointers, number, addNumber, max, tol,
                                       shift);
    }
  } else {
    this->computeEigenValuesUnsym(pointers, number, addNumber, max, tol,
                                     shift);
    
  }
}

void StaticSolutionState::computeEigenValuesSymUpper(
    PointerCollection &pointers, indexType number, indexType addNumber,
    bool max, prec tol, prec shift) {
#ifdef USE_SPECTRA

  auto &Logger = pointers.getSPDLogger();

  this->eigenValues.clear();
  this->eigenVectors.clear();

  if (addNumber > this->Kaa.cols())
    addNumber = static_cast<indexType>(this->Kaa.cols());

  Eigen::Matrix<prec, 1, Eigen::Dynamic> evalues;
  Eigen::Matrix<prec, Eigen::Dynamic, Eigen::Dynamic> evectors;

  if (max) {
    Spectra::SparseSymMatProd<prec, Eigen::Upper, Eigen::ColMajor, indexType>
        op(this->Kaa);
    Spectra::SymEigsSolver<Spectra::SparseSymMatProd<
        prec, Eigen::Upper, Eigen::ColMajor, indexType>>
        eigs(op, number, addNumber);
    eigs.init();
    try {
      // int nconv = eigs.compute(Spectra::SortRule::LargestMagn, 500, tol);
      eigs.compute(Spectra::SortRule::LargestMagn, 500, tol);
    } catch (const std::exception &e) {
      Logger.error("Eigenvalue solver failed: {}", e.what());
    }
    if (eigs.info() == Spectra::CompInfo::Successful) {
      evalues = eigs.eigenvalues();
      evectors = eigs.eigenvectors();
      for (auto i = 0; i < number; ++i) {
        this->eigenValues.push_back(evalues(i));
        if (this->eigenVectors.size() < i + 1)
          this->eigenVectors.emplace_back();
        for (auto j = 0; j < this->NumberOfActiveEquations; ++j) {
          this->eigenVectors[i].push_back(evectors(j, i));
        }
      }
    }
  } else {

    Spectra::SparseSymShiftSolve<prec, Eigen::Upper, Eigen::ColMajor, indexType>
        op(this->Kaa);
    Spectra::SymEigsShiftSolver<Spectra::SparseSymShiftSolve<
        prec, Eigen::Upper, Eigen::ColMajor, indexType>>
        eigs(op, number, addNumber, shift);

    eigs.init();
    try {
      // int nconv = eigs.compute(Spectra::SortRule::LargestMagn, 500, tol);
      eigs.compute(Spectra::SortRule::LargestMagn, 500, tol);
    } catch (const std::exception &e) {
      Logger.error("Eigenvalue solver failed: {}", e.what());
    }

    if (eigs.info() == Spectra::CompInfo::Successful) {
      evalues = eigs.eigenvalues();
      evectors = eigs.eigenvectors();
      for (auto i = 0; i < number; ++i) {
        this->eigenValues.push_back(evalues(i));
        if (this->eigenVectors.size() < i + 1)
          this->eigenVectors.emplace_back();
        for (auto j = 0; j < this->NumberOfActiveEquations; ++j) {
          this->eigenVectors[i].push_back(evectors(j, i));
        }
      }
    }
  }

  Logger.info(evalues);

#endif // USE_SPECTRA
}

void StaticSolutionState::computeEigenValuesSymLower(
    PointerCollection &pointers, indexType number, indexType addNumber,
    bool max, prec tol, prec shift)
{

#ifdef USE_SPECTRA
  auto &Logger = pointers.getSPDLogger();

  this->eigenValues.clear();
  this->eigenVectors.clear();

  if (addNumber > this->Kaa.cols())
    addNumber = static_cast<indexType>(this->Kaa.cols());

  Eigen::Matrix<prec, 1, Eigen::Dynamic> evalues;
  Eigen::Matrix<prec, Eigen::Dynamic, Eigen::Dynamic> evectors;

  if (max) {
    Spectra::SparseSymMatProd<prec, Eigen::Lower, Eigen::ColMajor, indexType>
        op(this->Kaa);
    Spectra::SymEigsSolver<Spectra::SparseSymMatProd<
        prec, Eigen::Lower, Eigen::ColMajor, indexType>>
        eigs(op, number, addNumber);
    eigs.init();
    try {
      // int nconv = eigs.compute(Spectra::SortRule::LargestMagn, 500, tol);
      eigs.compute(Spectra::SortRule::LargestMagn, 500, tol);
    } catch (const std::exception &e) {
      Logger.error("Eigenvalue solver failed: {}", e.what());
    }
    if (eigs.info() == Spectra::CompInfo::Successful) {
      evalues = eigs.eigenvalues();
      evectors = eigs.eigenvectors();
      for (auto i = 0; i < number; ++i) {
        this->eigenValues.push_back(evalues(i));
        if (this->eigenVectors.size() < i + 1)
          this->eigenVectors.emplace_back();
        for (auto j = 0; j < this->NumberOfActiveEquations; ++j) {
          this->eigenVectors[i].push_back(evectors(j, i));
        }
      }
    }
  } else {

    Spectra::SparseSymShiftSolve<prec, Eigen::Lower, Eigen::ColMajor, indexType>
        op(this->Kaa);
    Spectra::SymEigsShiftSolver<Spectra::SparseSymShiftSolve<
        prec, Eigen::Lower, Eigen::ColMajor, indexType>>
        eigs(op, number, addNumber, shift);

    eigs.init();
    try {
      // int nconv = eigs.compute(Spectra::SortRule::LargestMagn, 500, tol);
      eigs.compute(Spectra::SortRule::LargestMagn, 500, tol);
    } catch (const std::exception &e) {
      Logger.error("Eigenvalue solver failed: {}", e.what());
    }

    if (eigs.info() == Spectra::CompInfo::Successful) {
      evalues = eigs.eigenvalues();
      evectors = eigs.eigenvectors();
      for (auto i = 0; i < number; ++i) {
        this->eigenValues.push_back(evalues(i));
        if (this->eigenVectors.size() < i + 1)
          this->eigenVectors.emplace_back();
        for (auto j = 0; j < this->NumberOfActiveEquations; ++j) {
          this->eigenVectors[i].push_back(evectors(j, i));
        }
      }
    }
  }

  Logger.info(evalues);

#endif // USE_SPECTRA
}
void StaticSolutionState::computeEigenValuesUnsym(PointerCollection &pointers,
                                                  indexType number,
                                                  indexType addNumber, bool max,
                                                  prec tol, prec shift)
{

#ifdef USE_SPECTRA

  auto &Logger=pointers.getSPDLogger();

  this->eigenValues.clear();
  this->eigenVectors.clear();

  if (addNumber > this->Kaa.cols())
    addNumber = static_cast<indexType>(this->Kaa.cols());

  Eigen::Matrix<prec, 1, Eigen::Dynamic> evalues;
  Eigen::Matrix<prec, Eigen::Dynamic, Eigen::Dynamic> evectors;

  if (max) {
    Spectra::SparseGenMatProd<prec, Eigen::ColMajor, indexType>
        op(this->Kaa);
    Spectra::GenEigsSolver<Spectra::SparseGenMatProd<prec, Eigen::ColMajor, indexType>>
        eigs(op, number, addNumber);
    eigs.init();
    try {
      // int nconv = eigs.compute(Spectra::SortRule::LargestMagn, 500, tol);
      eigs.compute(Spectra::SortRule::LargestMagn, 500, tol);
    } catch (const std::exception &e) {
      Logger.error("Eigenvalue solver failed: {}", e.what());
    }
    if (eigs.info() == Spectra::CompInfo::Successful) {
      evalues = eigs.eigenvalues().real();
      evectors = eigs.eigenvectors().real();
      for (auto i = 0; i < number; ++i) {
        this->eigenValues.push_back(evalues(i));
        if (this->eigenVectors.size() < i + 1)
          this->eigenVectors.emplace_back();
        for (auto j = 0; j < this->NumberOfActiveEquations; ++j) {
          this->eigenVectors[i].push_back(evectors(j, i));
        }
      }
    }
  } else {

    Spectra::SparseGenRealShiftSolve<prec, Eigen::ColMajor, indexType>
        op(this->Kaa);
    Spectra::SymEigsShiftSolver<Spectra::SparseGenRealShiftSolve<
        prec, Eigen::ColMajor, indexType>>
        eigs(op, number, addNumber, shift);

    eigs.init();
    try {
      // int nconv = eigs.compute(Spectra::SortRule::LargestMagn, 500, tol);
      eigs.compute(Spectra::SortRule::LargestMagn, 500, tol);
    } catch (const std::exception &e) {
      Logger.error("Eigenvalue solver failed: {}", e.what());
    }

    if (eigs.info() == Spectra::CompInfo::Successful) {
      evalues = eigs.eigenvalues();
      evectors = eigs.eigenvectors();
      for (auto i = 0; i < number; ++i) {
        this->eigenValues.push_back(evalues(i));
        if (this->eigenVectors.size() < i + 1)
          this->eigenVectors.emplace_back();
        for (auto j = 0; j < this->NumberOfActiveEquations; ++j) {
          this->eigenVectors[i].push_back(evectors(j, i));
        }
      }
    }
  }

  Logger.info("Result of unsymmetric eigenvalue solver:");
  Logger.info(evalues);

#endif // USE_SPECTRA

};



auto StaticSolutionState::getSolution(GenericNodes &node)
    -> Types::Vector3<prec> {
  Types::Vector3<prec> sol;
  Types::Vector3<indexType> ids = {node.getDegreeOfFreedom(0).getId(),
                              node.getDegreeOfFreedom(1).getId(),
   node.getDegreeOfFreedom(2).getId()};
  sol = this->Solution(ids);
  return sol;
}

auto StaticSolutionState::getSolution(indexType globalId) -> prec {
  prec ret = this->Solution(globalId);

  return ret;
}

auto StaticSolutionState::getSolution(std::vector<DegreeOfFreedom *> &Dofs)
    -> Eigen::Matrix<prec, Eigen::Dynamic, 1> {

  auto Size = static_cast<indexType>(Dofs.size());
  Eigen::Matrix<prec, Eigen::Dynamic, 1> RetVec(Size);
  if (Size > 0) {
    for (auto i = 0; i < Size; ++i) {
      RetVec(i) = this->Solution(Dofs[i]->getId());
    }
  }

  return RetVec;
}

auto StaticSolutionState::getIncrementalSolution(
    std::vector<DegreeOfFreedom *> &Dofs)
    -> Eigen::Matrix<prec, Eigen::Dynamic, 1> {
  auto Size = static_cast<indexType>(Dofs.size());
  Types::VectorX<prec> RetVec(Size);
  if (Size > 0) {
    for (auto i = 0; i < Size; ++i) {
      RetVec(i) = this->IncSolution(Dofs[i]->getId());
    }
  }
  return RetVec;
}

auto StaticSolutionState::getNewtonSolution(
    std::vector<DegreeOfFreedom *> &Dofs)
    -> Eigen::Matrix<prec, Eigen::Dynamic, 1> {
  auto Size = static_cast<indexType>(Dofs.size());
  Types::VectorX<prec> RetVec(Size);
  if (Size > 0) {
    for (auto i = 0; i < Size; ++i) {
      RetVec(i) = this->NewtonSolution(Dofs[i]->getId());
    }
  }
  return RetVec;
}

void StaticSolutionState::resetSolution() {
  this->GenericSolutionState::resetSolution();
  this->Solution.setZero();
  this->IncSolution.setZero();
  this->dIncSolution.setZero();
  this->NewtonSolution.setZero();
}

void StaticSolutionState::printSpMat(PointerCollection &pointers) {
  auto infos = pointers.getInfoData();
  std::string directory = infos->fileNames[FileHandling::directory];
  std::string filename = infos->fileNames[FileHandling::infile];
  filename += ".ma";
  std::string temp = directory + filename;

  std::ofstream myfile;
  myfile.open(temp);
  myfile << std::setprecision(20);
  myfile << std::scientific;
  for (int k = 0; k < this->Kaa.outerSize(); ++k)
    for (typename Eigen::SparseMatrix<prec, 0, indexType>::InnerIterator it(
             this->Kaa, k);
         it; ++it) {
      myfile << it.row() + 1 << " " << it.col() + 1 << " " << it.value()
             << std::endl;
      it.value();
      it.row();   // row index
      it.col();   // col index (here it is equal to k)
      it.index(); // inner index, here it is equal to it.row()
    }

  myfile.close();
}

void StaticSolutionState::ctestout(PointerCollection &pointers) {

  auto &Logger = pointers.getSPDLogger();
  Logger.critical(this->Solution);

}
void StaticSolutionState::toFile(PointerCollection &pointers,
                                 std::ofstream &out) {
  GenericSolutionState::toFile(pointers, out);
  writeEigenMatrix(out, Solution);
  writeEigenMatrix(out, IncSolution);
  writeEigenMatrix(out, dIncSolution);
  writeEigenMatrix(out, NewtonSolution);
  writeEigenMatrix(out, Rhs);
  writeEigenMatrix(out, RhsB);
  writeEigenMatrix(out, eqSol);

  writeEigenSparseMatrix(out, Kaa);
  writeEigenSparseMatrix(out, Kab);
  writeEigenSparseMatrix(out, Kba);
  writeEigenSparseMatrix(out, Kbb);
}
void StaticSolutionState::fromFile(PointerCollection &pointers,
                                   std::ifstream &in) {
  GenericSolutionState::fromFile(pointers, in);
  readEigenMatrix(in, Solution);
  readEigenMatrix(in, IncSolution);
  readEigenMatrix(in, dIncSolution);
  readEigenMatrix(in, NewtonSolution);
  readEigenMatrix(in, Rhs);
  readEigenMatrix(in, RhsB);
  readEigenMatrix(in, eqSol);

  readEigenSparseMatrix(in, Kaa);
  readEigenSparseMatrix(in, Kab);
  readEigenSparseMatrix(in, Kba);
  readEigenSparseMatrix(in, Kbb);
  this->solver->analyze(Kaa);

  this->incSol.resize(this->NumberOfinActiveEquations);
  this->incSol.setZero();
}

void StaticSolutionState::RVEDatatoFile(PointerCollection &pointers,
                                        std::ofstream &out) {
  GenericSolutionState::RVEDatatoFile(pointers, out);
  writeEigenMatrix(out, Solution);
  writeEigenMatrix(out, IncSolution);
  writeEigenMatrix(out, dIncSolution);
  writeEigenMatrix(out, NewtonSolution);
  
}

void StaticSolutionState::RVEDatafromFile(PointerCollection &pointers,
                                          std::ifstream &in) {
  GenericSolutionState::RVEDatafromFile(pointers, in);
  readEigenMatrix(in, Solution);
  readEigenMatrix(in, IncSolution);
  readEigenMatrix(in, dIncSolution);
  readEigenMatrix(in, NewtonSolution);

  this->solver->analyze(Kaa);

  this->incSol.resize(this->NumberOfinActiveEquations);
  this->incSol.setZero();
}



void StaticSolutionState::computeConditionNumber(PointerCollection &pointers) {

  auto &Logger = pointers.getSPDLogger();

  auto numEq = this->Kaa.cols();
  Types::VectorX<prec> rr = Types::VectorX<prec>::Random(numEq);
  auto rrNorm = sqrt(rr.dot(rr));
  for (auto i = 0; i < 40; ++i) {
    rr /= rrNorm;
    rr = this->Kaa * rr;
    rrNorm = sqrt(rr.dot(rr));
  }
  auto evMax = rrNorm;

  rr = Eigen::Matrix<prec, Eigen::Dynamic, 1>::Random(numEq, 1);
  this->solver->solve(rr, rr);
  rrNorm = sqrt(rr.dot(rr));
  for (auto i = 0; i < 40; ++i) {
    rr /= rrNorm;
    this->solver->solve(rr, rr);
    rrNorm = sqrt(rr.dot(rr));
  }
  auto evMin = (prec)1 / rrNorm;

  prec cond = evMax / evMin;
  pointers.getCompuationData()->add("CPPFEMConditionNumber", cond);

  Logger.info("Result of Condition Number Computation:");
  Logger.info("    Minimum Eigenvalue:     {:>12.6e}",evMin);
  Logger.info("    Maximum Eigenvalue:     {:>12.6e}",evMax);
  Logger.info("    Condition Number:       {:>12.6e}",evMax/evMin);


  Eigen::SparseMatrix<prec, 0, indexType> TempMat, ScalMat;
  TempMat = this->Kaa;

  auto neqs = TempMat.rows();
  Eigen::Matrix<prec, Eigen::Dynamic, 1> tempVec(neqs);
  for (auto i = 0; i < neqs; ++i) {
    tempVec(i) = (prec)1 / sqrt(abs(TempMat.coeffRef(i, i)));
  }
  ScalMat = tempVec.asDiagonal();
  TempMat = ScalMat * TempMat * ScalMat;

  numEq = this->Kaa.cols();
  rr.resize(numEq);
  rr = Eigen::Matrix<prec, Eigen::Dynamic, 1>::Random(numEq, 1);
  rrNorm = sqrt(rr.dot(rr));
  for (auto i = 0; i < 20; ++i) {
    rr /= rrNorm;
    rr = TempMat * rr;
    rrNorm = sqrt(rr.dot(rr));
  }
  evMax = rrNorm;

  rr = Eigen::Matrix<prec, Eigen::Dynamic, 1>::Random(numEq, 1);
  this->solver->factorize(TempMat);
  this->solver->solve(rr, rr);
  rrNorm = sqrt(rr.dot(rr));
  for (auto i = 0; i < 20; ++i) {
    rr /= rrNorm;
    this->solver->solve(rr, rr);
    rrNorm = sqrt(rr.dot(rr));
  }
  evMin = (prec)1 / rrNorm;

  cond = evMax / evMin;
  pointers.getCompuationData()->add("CPPFEMConditionNumberScaled", cond);

  Logger.info("Result of Condition Number Computation:");
  Logger.info("    Minimum Eigenvalue:     {:>12.6e}",evMin);
  Logger.info("    Maximum Eigenvalue:     {:>12.6e}",evMax);
  Logger.info("    Condition Number:       {:>12.6e}",evMax/evMin);
  this->solver->factorize(this->Kaa);
}

} /* namespace HierAMuS */
