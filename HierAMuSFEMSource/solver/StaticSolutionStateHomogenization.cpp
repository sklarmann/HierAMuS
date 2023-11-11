// Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause


#include <memory>
#include <types/MatrixTypes.h>

#include <control/HandlingStructs.h>
#include <control/OutputHandler.h>
#include <pointercollection/pointercollection.h>

#include <solver/GenericSolutionState.h>
#include <solver/StaticSolutionStateHomogenization.h>

#include <control/HandlingStructs.h>
#include <pointercollection/pointercollection.h>

#include <finiteElements/ElementList.h>

#include <fstream>

#include "Homogenization/Homogenization2DSolid.h"
#include "solver/Homogenization/Homogenization2DSolid.h"
#include "solver/Homogenization/Homogenization3DSolid.h"
#include "solver/Homogenization/HomogenizationBeam.h"
#include "solver/Homogenization/HomogenizationShell.h"
#include "solver/Homogenization/Homogenization3DThermoMechBeam.h"

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

#include "control/BinaryWrite.h"


#include "spdlog/fmt/ostr.h"

#ifdef USE_SPECTRA
#include <Spectra/GenEigsRealShiftSolver.h>
#include <Spectra/GenEigsSolver.h>
#include <Spectra/MatOp/SparseGenMatProd.h>
#include <Spectra/MatOp/SparseGenRealShiftSolve.h>
#include <Spectra/MatOp/SparseSymMatProd.h>
#include <Spectra/MatOp/SparseSymShiftSolve.h>
#include <Spectra/SymEigsShiftSolver.h>
#include <Spectra/SymEigsSolver.h>
#endif // USE_SPECTRA

#include <iostream>

namespace HierAMuS {

StaticSolutionStateHomogenization::StaticSolutionStateHomogenization(
    ParameterList &parameter)
    : StaticSolutionState(parameter) {}

StaticSolutionStateHomogenization::StaticSolutionStateHomogenization(
    const StaticSolutionStateHomogenization &other) : StaticSolutionState(other)
{

  switch (other.homogenizationData->getType()) {
  case 1:
    this->homogenizationData = std::make_unique<Homogenization2DSolid>(reinterpret_cast<Homogenization2DSolid&>((*other.homogenizationData)));
    break;
  case 2:
    this->homogenizationData = std::make_unique<HomogenizationBeam>(
        reinterpret_cast<HomogenizationBeam &>((*other.homogenizationData)));
    break;
  case 3:
    this->homogenizationData = std::make_unique<Homogenization3DSolid>(
        reinterpret_cast<Homogenization3DSolid &>((*other.homogenizationData)));
    break;
  case 4:
    this->homogenizationData = std::make_unique<HomogenizationShell>(
        reinterpret_cast<HomogenizationShell &>((*other.homogenizationData)));
    break;
  case 20:
    this->homogenizationData = std::make_unique<Homogenization3DThermoMechBeam>(
        reinterpret_cast<Homogenization3DThermoMechBeam &>(
            (*other.homogenizationData)));
    break;
  }

  this->currStrains = other.currStrains;
  this->strainsIncrement = other.strainsIncrement;
  this->homogenizedData = other.homogenizedData;
}

StaticSolutionStateHomogenization::~StaticSolutionStateHomogenization() {
  this->Kaa.setZero();
}

auto StaticSolutionStateHomogenization::getCopy()
    -> std::shared_ptr<GenericSolutionState>  {


  std::shared_ptr<GenericSolutionState> tor =
      std::make_shared<StaticSolutionStateHomogenization>(*this);

  return tor;
}

void StaticSolutionStateHomogenization::setInitialValues(
    indexType numberOfEquations, indexType numberOfActiveEquations) {

  StaticSolutionState::setInitialValues(numberOfEquations,
                                        numberOfActiveEquations);
}

auto StaticSolutionStateHomogenization::residual() -> prec {
  return this->residualVal;
}

void StaticSolutionStateHomogenization::setEquationZero() {
  StaticSolutionState::setEquationZero();
}

void StaticSolutionStateHomogenization::solve(PointerCollection &pointers) {

  StaticSolutionState::solve(pointers);
}

void StaticSolutionStateHomogenization::computeLoads(
    PointerCollection &pointers) {


  this->homogenizationData->setPeriodicDisplacements(
      pointers, this->currStrains, this->strainsIncrement);

  auto strainDispInc = this->homogenizationData->getDisplacementIncrement(
      this->strainsIncrement);
  auto totdofs = pointers.getEquationHandler()->getNumberOfTotalEquations();
  for (auto i = 0; i < totdofs; ++i) {
    auto &Dof = pointers.getEquationHandler()->getDegreeOfFreedom(i);
    if (Dof.getStatus() == dofStatus::inactive) {
      // this->dIncSolution(Dof.getId()) += strainDispInc(Dof.getEqId());
      this->incSol(Dof.getEqId()) += strainDispInc(Dof.getEqId());
    }
  }
  StaticSolutionState::computeLoads(pointers);
}

void StaticSolutionStateHomogenization::assembleSystem(
    PointerCollection &pointers) {

  
  StaticSolutionState::assembleSystem(pointers);
  // ToDo Adding the strains
}

void StaticSolutionStateHomogenization::nextSolutionStep() {
  StaticSolutionState::nextSolutionStep();
}

void StaticSolutionStateHomogenization::updateSolution(
    PointerCollection &pointers) {
  StaticSolutionState::updateSolution(pointers);
  strainsIncrement.setZero();
}

void StaticSolutionStateHomogenization::resetSolution() {
  StaticSolutionState::resetSolution();
  this->currStrains.setZero();
  this->strainsIncrement.setZero();
}

void StaticSolutionStateHomogenization::setStrains(
    PointerCollection &pointers, Types::VectorX<prec> &strains) {
  if (this->currStrains.size() != strains.size()) {
    this->currStrains = strains;
    this->strainsIncrement = strains;
  } else {
    this->strainsIncrement = strains - this->currStrains;
    this->currStrains = strains;
  }

  auto &Logger = pointers.getSPDLogger();

  Logger.debug("  Updated strains to:     {}",this->currStrains.transpose());


}

auto StaticSolutionStateHomogenization::getCurrentStrains()
    -> Types::VectorX<prec> {
  return this->currStrains;
}

auto StaticSolutionStateHomogenization::getStrainsIncrement()
    -> Types::VectorX<prec> {
  return this->strainsIncrement;
}

void StaticSolutionStateHomogenization::initHomogenization(
    PointerCollection &pointers, indexType homogenizationType,
    ParameterList &parameters) {

  this->initHomogenizationType(pointers, homogenizationType);
  this->homogenizationData->init(pointers, parameters);
}

void StaticSolutionStateHomogenization::initHomogenizationType(
    PointerCollection &pointers, indexType homogenizationType) {
  switch (homogenizationType) {
  case 1:
    this->homogenizationData = std::make_unique<Homogenization2DSolid>();
    break;
  case 2:
    this->homogenizationData = std::make_unique<HomogenizationBeam>();
    break;
  case 3:
    this->homogenizationData = std::make_unique<Homogenization3DSolid>();
    break;
  case 4:
    this->homogenizationData = std::make_unique<HomogenizationShell>();
    break;
  case 20:
    this->homogenizationData =
        std::make_unique<Homogenization3DThermoMechBeam>();
    break;
  }
};

void StaticSolutionStateHomogenization::computeAMatrix(
    PointerCollection &pointers) {
  this->homogenizationData->computeAMatrix(pointers);
}

void StaticSolutionStateHomogenization::homogenize(
    PointerCollection &pointers) {
  auto &A = this->homogenizationData->getAMatrix();
  Types::MatrixXX<prec> temp = this->Kab * A;
  Types::MatrixXX<prec> temp2 = this->solver->solve(temp);
  Types::MatrixXX<prec> temp3 = A.transpose() * this->Kba;
  Types::MatrixXX<prec> tred = -temp3 * temp2;
  Types::MatrixXX<prec> tred2;
  if (this->symmetricSolver) {
    if (this->upper) {
      tred2 = A.transpose() * this->Kbb.selfadjointView<Eigen::Upper>() * A;
    } else {
      tred2 = A.transpose() * this->Kbb.selfadjointView<Eigen::Lower>() * A;
    }
  } else {
    tred2 = A.transpose() * this->Kbb * A;
  }

  

  tred += tred2;
  // A.transpose() * this->Kbb * A;
  tred /= this->homogenizationData->getDv();
  // sigmas
  Types::VectorX<prec> sig = -A.transpose() * RhsB - temp3 *eqSol;
  sig /= this->homogenizationData->getDv();

  prec threshold =
      tred.cwiseAbs().maxCoeff() * std::numeric_limits<prec>::epsilon() * 100;

  auto &Logger = pointers.getSPDLogger();
  Logger.debug("Material tangent:\n{}",tred.unaryExpr([threshold](prec d) {
                                    return abs(d) < threshold ? prec(0) : d;
                                  }));
  Logger.debug("Stresses:\n{}",sig.unaryExpr([threshold](prec d) {
                                    return abs(d) < threshold ? prec(0) : d;
                                  }));


  this->homogenizedData.C = tred;
  this->homogenizedData.sigma = sig;
}

auto StaticSolutionStateHomogenization::getHomogenizedData()
    -> GenericSolutionState::HomogenizedData {

  return this->homogenizedData;
}
void StaticSolutionStateHomogenization::toFile(PointerCollection &pointers,
                                               std::ofstream &out) {
  StaticSolutionState::toFile(pointers, out);
  writeEigenMatrix(out, this->currStrains);
  writeEigenMatrix(out, this->strainsIncrement);
  if (this->homogenizationData) {
    indexType homtype = this->homogenizationData->getType();
    writeScalar(out, homtype);
    this->homogenizationData->toFile(pointers, out);
  } else {
    indexType homtype = 0;
    writeScalar(out, homtype);
  }
}
void StaticSolutionStateHomogenization::fromFile(PointerCollection &pointers,
                                                 std::ifstream &in) {
  StaticSolutionState::fromFile(pointers, in);
  readEigenMatrix(in, this->currStrains);
  readEigenMatrix(in, this->strainsIncrement);
  indexType homtype;
  readScalar(in, homtype);
  if (homtype != 0) {
  }
  this->initHomogenizationType(pointers, homtype);
  this->homogenizationData->fromFile(pointers, in);
}

void StaticSolutionStateHomogenization::RVEDatatoFile(
    PointerCollection &pointers, std::ofstream &out) {
  StaticSolutionState::RVEDatatoFile(pointers, out);
  writeEigenMatrix(out, this->currStrains);
  writeEigenMatrix(out, this->strainsIncrement);
}

void StaticSolutionStateHomogenization::RVEDatafromFile(
    PointerCollection &pointers, std::ifstream &in) {
  StaticSolutionState::RVEDatafromFile(pointers, in);
  readEigenMatrix(in, this->currStrains);
  readEigenMatrix(in, this->strainsIncrement);
}

auto StaticSolutionStateHomogenization::getCMatrix() -> Types::MatrixXX<prec> {
  return this->homogenizedData.C;
}

auto StaticSolutionStateHomogenization::getStresses() -> Types::VectorX<prec> {
  return this->homogenizedData.sigma;
}

} /* namespace HierAMuS */
