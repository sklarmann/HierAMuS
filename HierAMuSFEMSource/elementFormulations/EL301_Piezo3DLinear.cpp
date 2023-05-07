// Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause



#include <elementFormulations/EL301_Piezo3DLinear.h>
#include <pointercollection/pointercollection.h>

#include <control/HandlingStructs.h>
#include <control/OutputHandler.h>

#include <equations/DegreeOfFreedom.h>
#include <equations/GenericNodes.h>
#include <equations/NodeSet.h>
#include <equations/Nodetypes.h>

#include <finiteElements/GenericFiniteElement.h>

#include <geometry/Base.h>


#include <materials/GenericMaterialFormulation.h>

#include "shapefunctions/IntegrationsPoints/dataClasses/GaussPoints.h"

#include <elementFormulations/GenericElementFormulation.h>
#include <solver/GenericSolutionState.h>

#include <types/MatrixTypes.h>

#include <Eigen/Eigenvalues>

namespace HierAMuS::Elementformulations {

EL301_Piezo3DLinear::EL301_Piezo3DLinear(PointerCollection *ptrCol)
    : GenericElementFormulation(ptrCol) {}

EL301_Piezo3DLinear::~EL301_Piezo3DLinear() = default;

void EL301_Piezo3DLinear::readData(PointerCollection &pointers,
                                   ParameterList &list) {
  this->meshIdDisp = list.getIndexVal("meshiddisp");
  this->meshIdPiezo = list.getIndexVal("meshidpiezo");

  this->intOrderDisp = list.getIndexVal("disporder");
  this->intOrderPiezo = list.getIndexVal("piezoorder");
  
  
  auto Log = pointers.getSPDLogger();

  Log.info("\n{:-<100}\n"
                "*   Element 301, specified Options\n"
                "    Mesh id for displacement nodes:         {:>12}\n"
                "    Shape function order for displacements: {:>12}\n"
                "    Mesh id for piezo nodes:                {:>12}\n"
                "    Shape function order for piezo:         {:>12}\n"
                "{:-<100}\n",
                "",
                this->meshIdDisp,
                this->intOrderDisp,
                this->meshIdPiezo,
                this->intOrderPiezo,
                "");


  this->messageUnprocessed(pointers, list,"EL301_Piezo3DLinear");
}

void EL301_Piezo3DLinear::setDegreesOfFreedom(
  PointerCollection& pointers, FiniteElement::GenericFiniteElement *elem) {

  elem->setH1Shapes(pointers, this->meshIdDisp, this->intOrderDisp);
  elem->setH1Shapes(pointers, this->meshIdPiezo, this->intOrderPiezo);
}

auto EL301_Piezo3DLinear::getDofs(PointerCollection& pointers, FiniteElement::GenericFiniteElement *elem) -> std::vector<DegreeOfFreedom*> {
  std::vector<DegreeOfFreedom *> Dofs, dispDofs, piezoDofs;
  elem->getH1Dofs(pointers, dispDofs, this->meshIdDisp,
                  this->intOrderDisp);
  elem->getH1Dofs(pointers, piezoDofs, this->meshIdPiezo,
                  this->intOrderPiezo);
  Dofs.insert(Dofs.begin(), dispDofs.begin(), dispDofs.end());
  Dofs.insert(Dofs.end(), piezoDofs.begin(), piezoDofs.end());
  return Dofs;
}

auto EL301_Piezo3DLinear::getNumberOfIntergrationPoints(
    PointerCollection &pointers, FiniteElement::GenericFiniteElement *elem)
    -> indexType {
  auto GP = elem->getIntegrationPoints(pointers);
  GP.setOrder(this->intOrderDisp * 2);
  return GP.getTotalGP();
}

void EL301_Piezo3DLinear::toParaviewAdaper(
    PointerCollection &pointers, FiniteElement::GenericFiniteElement *elem,
    vtkPlotInterface &paraviewAdapter, ParaviewSwitch control)
{

  pointers.getSPDLogger().warn("Plot functionality for EL301_Piezo3DLinear not implemented!");
}

void EL301_Piezo3DLinear::setTangentResidual(
  PointerCollection& pointers, FiniteElement::GenericFiniteElement *elem,
  Types::MatrixXX<prec> &stiffness, Types::VectorX<prec> &residual, std::vector<DegreeOfFreedom *> &Dofs) {

  Types::Matrix3X<prec> shapeDerivDisp, shapeDerivPiezo;
  Types::VectorX<prec> shapeDisp, shapePiezo;
  Types::MatrixXX<prec> Bmat;
  Types::MatrixXX<prec> jacobi;
  Types::MatrixXX<prec> material;
  Types::VectorX<prec> solution;

  Types::VectorX<prec> epsilon, sigma;
  epsilon.resize(9);
  sigma.resize(9);

  std::vector<DegreeOfFreedom *> dispDofs, piezoDofs;

  Dofs.clear();
  elem->getH1Dofs(pointers, dispDofs, this->meshIdDisp,
                  this->intOrderDisp);
  elem->getH1Dofs(pointers, piezoDofs, this->meshIdPiezo,
                  this->intOrderPiezo);
  Dofs.insert(Dofs.begin(), dispDofs.begin(), dispDofs.end());
  Dofs.insert(Dofs.end(), piezoDofs.begin(), piezoDofs.end());

  auto numDispDofs = static_cast<indexType>(dispDofs.size());
  auto numPiezoDofs = static_cast<indexType>(piezoDofs.size());

  auto numDofs = static_cast<indexType>(Dofs.size());

  stiffness.resize(numDofs, numDofs);
  stiffness.setZero();

  residual.resize(numDofs);
  residual.setZero();

  elem->getSolution(pointers, Dofs, solution);

  material.resize(9, 9);
  material.setZero();

  Bmat.resize(9, numDofs);
  Bmat.setZero();


  
  auto GP = elem->getIntegrationPoints(pointers);
  GP.setOrder(this->intOrderDisp * 2);

  prec da;

  Materials::MaterialTransferData materialData;

  for (auto i : GP) {
    jacobi = elem->getJacobian(pointers, i);
    auto H1ShapesDisp =
        elem->getH1Shapes(pointers, this->intOrderDisp, jacobi, i);
    auto H1ShapesPiezo =
        elem->getH1Shapes(pointers, this->intOrderPiezo, jacobi, i);

    for (auto l = 0; l < numDispDofs / 3; l++) {
      Bmat(0, l * 3) = H1ShapesDisp.shapeDeriv(0, l);
      Bmat(1, l * 3 + 1) = H1ShapesDisp.shapeDeriv(1, l);
      Bmat(2, l * 3 + 2) = H1ShapesDisp.shapeDeriv(2, l);

      Bmat(3, l * 3) = H1ShapesDisp.shapeDeriv(1, l);
      Bmat(3, l * 3 + 1) = H1ShapesDisp.shapeDeriv(0, l);

      Bmat(4, l * 3) = H1ShapesDisp.shapeDeriv(2, l);
      Bmat(4, l * 3 + 2) = H1ShapesDisp.shapeDeriv(0, l);

      Bmat(5, l * 3 + 1) = H1ShapesDisp.shapeDeriv(2, l);
      Bmat(5, l * 3 + 2) = H1ShapesDisp.shapeDeriv(1, l);
    }
    for (auto l = 0; l < numPiezoDofs / 3; l++) {
      Bmat(6, numDispDofs + l * 3) = H1ShapesPiezo.shapeDeriv(0, l);
      Bmat(7, numDispDofs + l * 3) = H1ShapesPiezo.shapeDeriv(1, l);
      Bmat(8, numDispDofs + l * 3) = H1ShapesPiezo.shapeDeriv(2, l);
    }
    materialData.strains = Bmat * solution;
    elem->getMaterialFormulation(pointers)->getMaterialData(pointers, materialData,i);
    da = i.weight * jacobi.determinant();
    stiffness += Bmat.transpose() * materialData.materialTangent * Bmat * da;
    residual += Bmat.transpose() * materialData.stresses * da;
  }

  // for (auto i = 0; i < xsi.size(); ++i) {
  //   elem->getJacobian(*this->ptrCol, jacobi, xsi[i], eta[i]);
  //   elem->getH1Shapes(*(this->ptrCol), this->intOrderDisp, jacobi, shapeT,
  //                     shapeDeriv, xsi[i], eta[i]);
  //   for (auto j = 0; j < numDofs / 3; ++j) {
  //     Bmat(0, j * 3) = shapeDeriv(0, j);
  //     Bmat(1, j * 3 + 1) = shapeDeriv(1, j);
  //     Bmat(2, j * 3) = shapeDeriv(1, j);
  //     Bmat(2, j * 3 + 1) = shapeDeriv(0, j);
  //     Bmat(3, j * 3 + 2) = shapeDeriv(0, j);
  //     Bmat(4, j * 3 + 2) = shapeDeriv(1, j);
  //   }

  //   da = (jacobi.determinant()) * weight[i];
  //   stiffness += Bmat.transpose() * material * Bmat * da;
  // }
  // residual = stiffness * solution;
}

  /**
void EL301_Piezo3DLinear::projectedToParaview(
    PointerCollection &pointers, FiniteElement::GenericFiniteElement *elem,
    vtkPlotInterface &paraviewAdapter) {

  Types::Matrix3X<prec> shapeDerivDisp, shapeDerivPiezo, shapeDerivInterp;
  Types::VectorX<prec> shapeDisp, shapePiezo, shapeInterp;
  Types::MatrixXX<prec> jacobi;
  Types::MatrixXX<prec> material;
  Types::VectorX<prec> solution;

  Types::VectorX<prec> epsilon, sigma, gpDetJDummy;
  gpDetJDummy.resize(1);
  epsilon.resize(9);
  sigma.resize(9);

  gpDetJDummy(0) = 1;

  std::vector<DegreeOfFreedom *> dispDofs, piezoDofs, Dofs;
  elem->getH1Dofs(pointers, dispDofs, this->meshIdDisp,
                  this->intOrderDisp);
  elem->getH1Dofs(pointers, piezoDofs, this->meshIdPiezo,
                  this->intOrderPiezo);
  Dofs.insert(Dofs.begin(), dispDofs.begin(), dispDofs.end());
  Dofs.insert(Dofs.end(), piezoDofs.begin(), piezoDofs.end());

  auto numDispDofs = static_cast<indexType>(dispDofs.size());
  auto numPiezoDofs = static_cast<indexType>(piezoDofs.size());
  auto numDofs = static_cast<indexType>(Dofs.size());

  elem->getSolution(pointers, Dofs, solution);

  Types::MatrixXX<prec> Bmat;
  Bmat.resize(9, numDofs);
  Bmat.setZero();
  std::vector<prec> xsi, eta, zeta, weight;
  elem->getGaussPoints(this->intOrderDisp + 1, weight, xsi, eta, zeta);
  auto GP = elem->getIntegrationPoints(pointers);
  GP.setOrder(this->intOrderDisp * 2);

  Materials::MaterialTransferData materialData;

  for (auto i: GP) {
    jacobi = elem->getJacobian(pointers, i);
    auto H1ShapesDisp =
        elem->getH1Shapes(pointers, this->intOrderDisp, jacobi, i);
    auto H1ShapesPiezo =
        elem->getH1Shapes(pointers, this->intOrderPiezo, jacobi, i);
    

    for (auto l = 0; l < numDispDofs / 3; l++) {
      Bmat(0, l * 3) = H1ShapesDisp.shapeDeriv(0, l);
      Bmat(1, l * 3 + 1) = H1ShapesDisp.shapeDeriv(1, l);
      Bmat(2, l * 3 + 2) = H1ShapesDisp.shapeDeriv(2, l);

      Bmat(3, l * 3) = H1ShapesDisp.shapeDeriv(1, l);
      Bmat(3, l * 3 + 1) = H1ShapesDisp.shapeDeriv(0, l);

      Bmat(4, l * 3) = H1ShapesDisp.shapeDeriv(2, l);
      Bmat(4, l * 3 + 2) = H1ShapesDisp.shapeDeriv(0, l);

      Bmat(5, l * 3 + 1) = H1ShapesDisp.shapeDeriv(2, l);
      Bmat(5, l * 3 + 2) = H1ShapesDisp.shapeDeriv(1, l);
    }
    for (auto l = 0; l < numPiezoDofs / 3; l++) {
      Bmat(6, numDispDofs + l * 3) = H1ShapesPiezo.shapeDeriv(0, l);
      Bmat(7, numDispDofs + l * 3) = H1ShapesPiezo.shapeDeriv(1, l);
      Bmat(8, numDispDofs + l * 3) = H1ShapesPiezo.shapeDeriv(2, l);
    }
    materialData.strains = Bmat * solution;
    elem->getMaterialFormulation(pointers)->getMaterialData(pointers, materialData,i);
    elem->projectOnVertsParaview(pointers, paraviewAdapter, gpDetJDummy, xsi[i],
                                 eta[i], zeta[i], weight[i],
                                 paraviewNames::weightName());
    elem->projectOnVertsParaview(pointers, paraviewAdapter,
                                 materialData.strains, xsi[i], eta[i], zeta[i],
                                 weight[i], paraviewNames::strainName());
    elem->projectOnVertsParaview(pointers, paraviewAdapter,
                                 materialData.stresses, xsi[i], eta[i], zeta[i],
                                 weight[i], paraviewNames::stressName());
  }
}
**/
} // namespace HierAMuS

