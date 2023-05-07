// Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause





#include <elementFormulations/EL202_Piezo2D.h>

#include <control/HandlingStructs.h>
#include <control/OutputHandler.h>
#include <pointercollection/pointercollection.h>

#include <equations/DegreeOfFreedom.h>
#include <equations/GenericNodes.h>
#include <equations/NodeSet.h>
#include <equations/Nodetypes.h>

#include <finiteElements/GenericFiniteElement.h>

#include <geometry/Base.h>


#include "shapefunctions/IntegrationsPoints/dataClasses/GaussPoints.h"

#include <elementFormulations/GenericElementFormulation.h>
#include <solver/GenericSolutionState.h>

#include <types/MatrixTypes.h>

#include <Eigen/Eigenvalues>

namespace HierAMuS::Elementformulations {

EL202_Piezo2D::EL202_Piezo2D(PointerCollection *ptrCol)
    : GenericElementFormulation(ptrCol) {}

EL202_Piezo2D::~EL202_Piezo2D() = default;

void EL202_Piezo2D::readData(PointerCollection &pointers, ParameterList &list) {

  this->meshIdDisp = list.getIndexVal("meshiddisp");
  this->intOrderDisp = list.getIndexVal("disporder");

  this->emod = list.getPrecVal("emodul");
  this->nu = list.getPrecVal("nu");

  this->mu1 = list.getPrecVal("mu1");
  this->mu2 = list.getPrecVal("mu2");

  this->e31 = list.getPrecVal("e31");
  this->e33 = list.getPrecVal("e33");
  this->e15 = list.getPrecVal("e15");


  auto Log = pointers.getSPDLogger();

  Log.info("\n{:-<100}\n"
                "*   Element 201, specified Options\n"
                "    Mesh id for displacement nodes:         {:>12}\n"
                "    Shape function order for displacements: {:>12}\n"
                "{:-<100}\n",
                "",
                this->meshIdDisp,
                this->intOrderDisp,
                "");

  this->messageUnprocessed(pointers, list, "EL202_Piezo2D");
}

void EL202_Piezo2D::setDegreesOfFreedom(
  PointerCollection& pointers, FiniteElement::GenericFiniteElement *elem) {

  elem->setH1Shapes(pointers, this->meshIdDisp, this->intOrderDisp);
}

auto EL202_Piezo2D::getDofs(PointerCollection& pointers, FiniteElement::GenericFiniteElement *elem) -> std::vector<DegreeOfFreedom *>
{
  std::vector<DegreeOfFreedom *> Dofs;
  elem->getH1Dofs(pointers, Dofs, this->meshIdDisp, this->intOrderDisp);
  return Dofs;
}

auto EL202_Piezo2D::getNumberOfIntergrationPoints(
    PointerCollection &pointers, FiniteElement::GenericFiniteElement *elem)
    -> indexType {

  auto GP = elem->getIntegrationPoints(pointers);
  GP.setOrder(this->intOrderDisp * 2);
  return GP.getTotalGP();
}

void EL202_Piezo2D::toParaviewAdaper(PointerCollection &pointers,
                                     FiniteElement::GenericFiniteElement *elem,
                                     vtkPlotInterface &paraviewAdapter,
                                     ParaviewSwitch control)
{
  pointers.getSPDLogger().warn("Plot functionality for element formulation EL202_Piezo2D not implemented!");
}

void EL202_Piezo2D::setTangentResidual(
  PointerCollection& pointers, FiniteElement::GenericFiniteElement *elem,
  Types::MatrixXX<prec> &stiffness, Types::VectorX<prec> &residual, std::vector<DegreeOfFreedom *> &Dofs) {
  
  Types::MatrixXX<prec> Bmat;
  Types::MatrixXX<prec> material;
  Types::VectorX<prec> solution;

  Dofs.clear();
  elem->getH1Dofs(pointers, Dofs, this->meshIdDisp, this->intOrderDisp);

  auto numDofs = static_cast<indexType>(Dofs.size());

  stiffness.resize(numDofs, numDofs);
  stiffness.setZero();

  residual.resize(numDofs);
  residual.setZero();

  elem->getSolution(pointers, Dofs, solution);

  material.resize(5, 5);
  material.setZero();
  prec fac = this->emod / ((prec)1 + this->nu) / ((prec)1 - (prec)2 * this->nu);
  material(0, 0) = ((prec)1 - this->nu);
  material(0, 1) = (this->nu);
  material(1, 1) = ((prec)1 - this->nu);
  material(1, 0) = (this->nu);
  material *= fac;
  material(2, 2) = this->emod / ((prec)1 + this->nu) / (prec)2;

  material(3, 3) = this->mu1;
  material(4, 4) = this->mu2;

  material(4, 0) = -this->e31;
  material(4, 1) = -this->e33;
  material(3, 2) = -this->e15;

  material(0, 4) = -this->e31;
  material(1, 4) = -this->e33;
  material(2, 3) = -this->e15;


  

  Bmat.resize(5, numDofs);
  Bmat.setZero();
  
  auto GP = elem->getIntegrationPoints(pointers);
  GP.setOrder(this->intOrderDisp * 2);

  for (auto ip : GP)
  {
    auto jaco = elem->getJacobian(pointers, ip);
    auto shapes = elem->getH1Shapes(pointers, this->intOrderDisp, jaco, ip);
    for (auto j = 0; j < numDofs / 3; ++j) {
      Bmat(0, j * 3)     = shapes.shapeDeriv(0, j);
      Bmat(1, j * 3 + 1) = shapes.shapeDeriv(1, j);
      Bmat(2, j * 3)     = shapes.shapeDeriv(1, j);
      Bmat(2, j * 3 + 1) = shapes.shapeDeriv(0, j);
      Bmat(3, j * 3 + 2) = shapes.shapeDeriv(0, j);
      Bmat(4, j * 3 + 2) = shapes.shapeDeriv(1, j);
    }

    
    auto da = jaco.determinant() * ip.weight;
    stiffness += Bmat.transpose() * material * Bmat * da;
  }
  
  residual = stiffness * solution;
}

} // namespace HierAMuS

