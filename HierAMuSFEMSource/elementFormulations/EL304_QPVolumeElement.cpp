// Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause




#include <elementFormulations/EL304_QPVolumeElement.h>
#include <pointercollection/pointercollection.h>

#include <control/HandlingStructs.h>
#include <control/OutputHandler.h>
#include "control/ParameterList.h"


#include <finiteElements/Volume.h>



#include <materials/GenericMaterialFormulation.h>

#include "shapefunctions/IntegrationsPoints/dataClasses/GaussPoints.h"

#include <elementFormulations/GenericElementFormulation.h>
#include <solver/GenericSolutionState.h>

#include <types/MatrixTypes.h>

#include <Eigen/Eigenvalues>

namespace HierAMuS {
namespace Elementformulations {

EL304_QPVolumeElement::EL304_QPVolumeElement(PointerCollection *ptrCol)
    : GenericElementFormulationInterface(ptrCol) {}

EL304_QPVolumeElement::~EL304_QPVolumeElement() {}

void EL304_QPVolumeElement::readData(PointerCollection &pointers,
                                           ParameterList &list) {
  this->dispOrder = list.getIndexVal("disporder");
  this->pressureOrder = this->dispOrder-1;
  this->meshIdDisp = list.getIndexVal("meshiddisp");

  this->mu = list.getPrecVal("mu");
  this->kappa = list.getPrecVal("kappa");
  

  auto &Log = pointers.getSPDLogger();

  Log.info("\n{:-<100}\n"
                "*   Element 304, specified Options\n"
                "    Mesh id for displacement nodes:         {:>12}\n"
                "    Shape function order for displacements: {:>12}\n"
                "    Shape function order for pressure:      {:>12}\n"
                "    Material parameters:\n"
                "      mu:                                   {:>12.4e}\n"
                "      kappa:                                {:>12.4e}\n"
                "{:-<100}\n",
                "",
                this->meshIdDisp,
                this->dispOrder,
                this->pressureOrder,
                this->mu,
                this->kappa,
                "");


  this->messageUnprocessed(pointers, list, "EL303_ThermoMechanikSolid3D");
}

void EL304_QPVolumeElement::setDegreesOfFreedom(
    PointerCollection &pointers, FiniteElement::Volume &elem) {

  elem.setH1Shapes(pointers, this->meshIdDisp, this->dispOrder);
}

void EL304_QPVolumeElement::AdditionalOperations(
  PointerCollection& pointers, FiniteElement::Volume &elem) {

}

auto EL304_QPVolumeElement::getDofs(PointerCollection& pointers, FiniteElement::GenericFiniteElement *elem) -> std::vector<DegreeOfFreedom*>{
  std::vector<DegreeOfFreedom *> ret, temp;

  elem->getH1Dofs(pointers, ret, this->meshIdDisp, this->dispOrder);
  return ret;
}

void EL304_QPVolumeElement::setTangentResidual(
  PointerCollection& pointers, FiniteElement::Volume &elem,
  Types::MatrixXX<prec> &stiffness, Types::VectorX<prec> &residual, std::vector<DegreeOfFreedom *> &Dofs) {

  this->setTangentResidualNonLinear(pointers, elem, stiffness, residual, Dofs);

}


void EL304_QPVolumeElement::setTangentResidualNonLinear(
  PointerCollection& pointers,
  FiniteElement::Volume &elem,
  Eigen::Matrix<prec, Eigen::Dynamic, Eigen::Dynamic> &stiffness,
  Eigen::Matrix<prec, Eigen::Dynamic, 1> &residual, std::vector<DegreeOfFreedom *> &Dofs) {

  indexType pressureDofs;
  if (pressureOrder == 0) {
    pressureDofs = 1;
  } else if(pressureOrder == 1){
    pressureDofs = 4;
  } else {
    throw std::runtime_error("Pressure order not supported!");
  }

  Dofs.clear();
  elem.getH1Dofs(pointers, Dofs, meshIdDisp, dispOrder);
  auto GP = this->getIntegrationPoints(pointers, elem);

  Types::Matrix66<prec> xIP;
  xIP.setZero();


  constexpr prec twoThird = prec(2) / prec(3);
  constexpr prec oneThird = prec(1) / prec(3);
  constexpr prec oneHalf = prec(1) / prec(2);

  xIP(0, 0) = twoThird;
  xIP(0, 1) = -oneThird;
  xIP(0, 2) = -oneThird;
  xIP(1, 0) = -oneThird;
  xIP(1, 1) = twoThird;
  xIP(1, 2) = -oneThird;
  xIP(2, 0) = -oneThird;
  xIP(2, 1) = -oneThird;
  xIP(2, 2) = twoThird;
  xIP(3, 3) = oneHalf;
  xIP(4, 4) = oneHalf;
  xIP(5, 5) = oneHalf;

  Types::MatrixXX<prec> xKm(Dofs.size(), Dofs.size()), xLm(Dofs.size(),pressureDofs), xC(pressureDofs,pressureDofs);
  Types::VectorX<prec> vf_in(Dofs.size());
  Types::VectorXT<prec> xNp(pressureDofs);
  xNp(0) = prec(1);
  xKm.setZero();
  vf_in.setZero();
  xLm.setZero();
  xC.setZero();

  Types::Matrix6X<prec> Bd = Types::Matrix6X<prec>::Zero(6, Dofs.size());
  Types::VectorXT<prec> Bv = Types::VectorXT<prec>::Zero(Dofs.size());

  stiffness.resize(Dofs.size(), Dofs.size());
  residual.resize(Dofs.size());
  stiffness.setZero();
  residual.setZero();

  Types::VectorX<prec> solution = elem.getSolution(pointers, Dofs);

  

  for (auto &iGP : GP) {
    Types::Matrix33<prec> jaco =  elem.getJacobian(pointers, iGP);
    auto shapes = elem.getH1Shapes(pointers, dispOrder, jaco, iGP);
    this->getBdBv(pointers, Bd, Bv, shapes);
    if (pressureDofs > 1) {
      xNp(1) = iGP.xi;
      xNp(2) = iGP.eta;
      xNp(3) = iGP.zeta;
    }
    prec detJ = jaco.determinant();
    prec dv = detJ * iGP.weight;
    prec dv1 = dv * prec(2) * this->mu;
    xKm += Bd.transpose() * xIP * Bd * dv1;
    xLm += Bv.transpose() * xNp * dv;
    dv /= this->kappa;
    xC += xNp.transpose() * xNp * dv;
  }
  Types::MatrixXX<prec> xCinv = xC.inverse();
  auto vp = xCinv * xLm.transpose() * solution;

  for (auto &iGP : GP) {
    Types::Matrix33<prec> jacobi = elem.getJacobian(pointers, iGP);
    auto shapes = elem.getH1Shapes(pointers, dispOrder, jacobi, iGP);
    prec detJ = jacobi.determinant();
    prec dv = detJ * iGP.weight;
    if (pressureDofs > 1) {
      xNp(1) = iGP.xi;
      xNp(2) = iGP.eta;
      xNp(3) = iGP.zeta;
    }
    prec pressure = xNp * vp;
    Types::Vector6<prec> eps = Bd * solution;
    Types::Vector6<prec> sig = xIP * eps * prec(2) * mu;
    residual += Bd.transpose() * sig * dv;
    residual += Bv.transpose() * pressure * dv;
  }
  stiffness = xKm + xLm * xCinv * xLm.transpose();
  
}



auto EL304_QPVolumeElement::getIntegrationPoints(
    PointerCollection &pointers, FiniteElement::GenericFiniteElement &elem)
    -> IntegrationPoints {
  auto GP = elem.getIntegrationPoints(pointers);
  GP.setOrder(dispOrder * 2);
  return GP;
}

void EL304_QPVolumeElement::getBdBv(PointerCollection &pointers,
                                    Types::Matrix6X<prec> &Bd,
                                    Types::VectorXT<prec> &Bv,
                                    Geometry::H1Shapes &shapes) {


  for (auto ii = 0; ii < shapes.shapes.rows(); ++ii) {
    Bd(0, ii*3 + 0) = shapes.shapeDeriv(0, ii);
    Bd(1, ii*3 + 1) = shapes.shapeDeriv(1, ii);
    Bd(2, ii*3 + 2) = shapes.shapeDeriv(2, ii);
    Bd(3, ii*3 + 0) = shapes.shapeDeriv(1, ii);
    Bd(3, ii*3 + 1) = shapes.shapeDeriv(0, ii);
    Bd(4, ii*3 + 0) = shapes.shapeDeriv(2, ii);
    Bd(4, ii*3 + 2) = shapes.shapeDeriv(0, ii);
    Bd(5, ii*3 + 1) = shapes.shapeDeriv(2, ii);
    Bd(5, ii*3 + 2) = shapes.shapeDeriv(1, ii);

    Bv(ii * 3) = shapes.shapeDeriv(0, ii);
    Bv(ii * 3 + 1) = shapes.shapeDeriv(1, ii);
    Bv(ii * 3 + 2) = shapes.shapeDeriv(2, ii);
  }

  
}

void EL304_QPVolumeElement::toParaviewAdaper(
    PointerCollection &pointers, FiniteElement::Volume &elem,
    vtkPlotInterface &paraviewAdapter, ParaviewSwitch control) {
  int matNum = static_cast<int>(elem.getMaterial()->getNumber());
  switch (control) {
  case ParaviewSwitch::Mesh: {
    elem.geometryToParaview(pointers, paraviewAdapter, 0, matNum);

  } break;
  case ParaviewSwitch::Solution: {
    elem.H1SolutionToParaview(pointers, paraviewAdapter, 0, matNum,
                               this->meshIdDisp, this->dispOrder,
                               paraviewNames::DisplacementName());
  } break;
  case ParaviewSwitch::Weights: {
    elem.computeWeightsParaview(pointers, paraviewAdapter, 0, matNum);
  } break;
  case ParaviewSwitch::ProjectedValues: {
  }break;
  default:
    break;
  }
}

auto EL304_QPVolumeElement::getHistoryDataStructure()
    -> const HistoryDataStructure & {
  return m_HistoryDataStructure;
}

auto EL304_QPVolumeElement::getNumberOfIntergrationPoints(
  PointerCollection& pointers, FiniteElement::GenericFiniteElement *elem) -> indexType {
  auto GP = this->getIntegrationPoints(pointers, *elem);
  return GP.getTotalGP();
}

const HistoryDataStructure EL304_QPVolumeElement::m_HistoryDataStructure({},{});

} // namespace Elementformulations
} // namespace HierAMuS
