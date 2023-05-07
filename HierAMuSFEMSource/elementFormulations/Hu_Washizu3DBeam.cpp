// Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause




#include <elementFormulations/Hu_Washizu3DBeam.h>
#include <pointercollection/pointercollection.h>

#include <control/HandlingStructs.h>
#include <control/OutputHandler.h>
#include <control/stringCommandHandler.h>

#include <equations/DegreeOfFreedom.h>
#include <equations/GenericNodes.h>
#include <equations/NodeSet.h>
#include <equations/Nodetypes.h>

#include <finiteElements/GenericFiniteElement.h>

#include <geometry/Base.h>

#include <math/Userconstants.h>

#include <materials/GenericMaterialFormulation.h>

#include "shapefunctions/IntegrationsPoints/GaussPoints.h"

#include <elementFormulations/GenericElementFormulation.h>
#include <solver/GenericSolutionState.h>

#include <MatrixTypes.h>

#include <Eigen/Eigenvalues>

namespace HierAMuS::Elementformulations {

template <typename prec, typename indexType>
Hu_Washizu3DBeam<prec, indexType>::Hu_Washizu3DBeam(
    PointerCollection<prec, indexType> *ptrCol)
    : GenericElementFormulation<prec, indexType>(ptrCol) {}

template <typename prec, typename indexType>
Hu_Washizu3DBeam<prec, indexType>::~Hu_Washizu3DBeam() {}

template <typename prec, typename indexType>
void Hu_Washizu3DBeam<prec, indexType>::readData(stringCommandHandler &Command,
                                                 Userconstants<prec> *ucons) {
  std::string temp;
  temp = Command.getRhs("meshiddisp");
  this->meshIdDisp = static_cast<indexType>(ucons->process(temp));
  temp = Command.getRhs("meshidrot");
  this->meshIdRot = static_cast<indexType>(ucons->process(temp));

  temp = Command.getRhs("disporder");
  this->intOrderDisp = static_cast<indexType>(ucons->process(temp));
  temp = Command.getRhs("rotorder");
  this->intOrderRot = static_cast<indexType>(ucons->process(temp));

  OutputHandler *Log;
  Log = &this->ptrCol->getInfoData()->Log;
  (*Log)(LogLevel::BasicLog, LogLevel::BasicLog)
      << "Hu_Washizu3DBeam, specified Options" << std::endl
      << "Mesh id for displacement nodes: " << this->meshIdDisp << std::endl
      << "Shape function order for displacements: " << this->intOrderDisp
      << "Mesh id for Rotation nodes: " << this->meshIdRot << std::endl
      << "Shape function order for rotation:" << this->intOrderRot << std::endl;
}

template <typename prec, typename indexType>
void Hu_Washizu3DBeam<prec, indexType>::setDegreesOfFreedom(
    FiniteElement::GenericFiniteElement<prec, indexType> *elem) {

  elem->setH1Shapes(*this->ptrCol, this->meshIdDisp, this->intOrderDisp);
  elem->setH1Shapes(*this->ptrCol, this->meshIdRot, this->intOrderRot);
}

template <typename prec, typename indexType>
void Hu_Washizu3DBeam<prec, indexType>::setTangentResidual(
    FiniteElement::GenericFiniteElement<prec, indexType> *elem, Types::MatrixXX<prec> &stiffness,
    Types::VectorX<prec> &residual,
    std::vector<DegreeOfFreedom<prec, indexType> *> &Dofs) {

  Types::Matrix3X<prec> shapeDerivDisp, shapeDerivPiezo;
  Types::VectorX<prec> shapeDisp, shapePiezo;
  Types::MatrixXX<prec> Bmat;
  Types::Matrix33<prec> jacobi;
  Types::MatrixXX<prec> material;
  Types::VectorX<prec> solution;

  Types::VectorX<prec> epsilon, sigma;
  epsilon.resize(6);
  sigma.resize(6);

  std::vector<DegreeOfFreedom<prec, indexType> *> dispDofs, piezoDofs;

  Dofs.clear();
  elem->getH1Dofs(*this->ptrCol, dispDofs, this->meshIdDisp,
                  this->intOrderDisp);
  elem->getH1Dofs(*this->ptrCol, piezoDofs, this->meshIdPiezo,
                  this->intOrderPiezo);
  Dofs.insert(Dofs.begin(), dispDofs.begin(), dispDofs.end());
  Dofs.insert(Dofs.end(), piezoDofs.begin(), piezoDofs.end());

  indexType numDispDofs = static_cast<indexType>(dispDofs.size());
  indexType numPiezoDofs = static_cast<indexType>(piezoDofs.size());

  indexType numDofs = static_cast<indexType>(Dofs.size());

  stiffness.resize(numDofs, numDofs);
  stiffness.setZero();

  residual.resize(numDofs);
  residual.setZero();

  elem->getSolution(*this->ptrCol, Dofs, solution);

  material.resize(9, 9);
  material.setZero();

  Bmat.resize(9, numDofs);
  Bmat.setZero();

  indexType approxGP = (this->intOrderDisp + 1) * (this->intOrderDisp + 1) *
                  (this->intOrderDisp + 1);

  std::vector<prec> xsi, eta, zeta, weight;
  xsi.reserve(approxGP), eta.reserve(approxGP), zeta.reserve(approxGP),
      weight.reserve(approxGP);
  elem->getGaussPoints(this->intOrderDisp + 1, weight, xsi, eta, zeta);
  prec da;
  for (auto i = 0; i < xsi.size(); i++) {
    elem->getJacobian(*this->ptrCol, jacobi, xsi[i], eta[i], zeta[i]);
    elem->getH1Shapes(*this->ptrCol, this->intOrderDisp, jacobi, shapeDisp,
                      shapeDerivDisp, xsi[i], eta[i], zeta[i]);

    elem->getH1Shapes(*this->ptrCol, this->intOrderPiezo, jacobi, shapePiezo,
                      shapeDerivPiezo, xsi[i], eta[i], zeta[i]);

    for (auto l = 0; l < numDispDofs / 3; l++) {
      Bmat(0, l * 3) = shapeDerivDisp(0, l);
      Bmat(1, l * 3 + 1) = shapeDerivDisp(1, l);
      Bmat(2, l * 3 + 2) = shapeDerivDisp(2, l);

      Bmat(3, l * 3) = shapeDerivDisp(1, l);
      Bmat(3, l * 3 + 1) = shapeDerivDisp(0, l);

      Bmat(4, l * 3) = shapeDerivDisp(2, l);
      Bmat(4, l * 3 + 2) = shapeDerivDisp(0, l);

      Bmat(5, l * 3 + 1) = shapeDerivDisp(2, l);
      Bmat(5, l * 3 + 2) = shapeDerivDisp(1, l);
    }
    for (auto l = 0; l < numPiezoDofs / 3; l++) {
      Bmat(6, numDispDofs + l * 3) = shapeDerivPiezo(0, l);
      Bmat(7, numDispDofs + l * 3) = shapeDerivPiezo(1, l);
      Bmat(8, numDispDofs + l * 3) = shapeDerivPiezo(2, l);
    }
    epsilon = Bmat * solution;
    elem->getMaterialFormulation()->getMaterialData(epsilon, sigma, material);
    da = weight[i] * jacobi.determinant();
    stiffness += Bmat.transpose() * material * Bmat * da;
    residual += Bmat.transpose() * sigma * da;
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

template <typename prec, typename indexType>
void Hu_Washizu3DBeam<prec, indexType>::projectedToParaview(
    PointerCollection<prec, indexType> &pointers,
    FiniteElement::GenericFiniteElement<prec, indexType> *elem, managementClass &paraviewAdapter) {

  Types::Matrix3X<prec> shapeDerivDisp, shapeDerivPiezo, shapeDerivInterp;
  Types::VectorX<prec> shapeDisp, shapePiezo, shapeInterp;
  Types::Matrix33<prec> jacobi;
  Types::MatrixXX<prec> material;
  Types::VectorX<prec> solution;

  Types::VectorX<prec> epsilon, sigma, gpDetJDummy;
  gpDetJDummy.resize(1);
  epsilon.resize(9);
  sigma.resize(9);

  gpDetJDummy(0) = 1;

  std::vector<DegreeOfFreedom<prec, indexType> *> dispDofs, piezoDofs, Dofs;
  elem->getH1Dofs(*this->ptrCol, dispDofs, this->meshIdDisp,
                  this->intOrderDisp);
  elem->getH1Dofs(*this->ptrCol, piezoDofs, this->meshIdPiezo,
                  this->intOrderPiezo);
  Dofs.insert(Dofs.begin(), dispDofs.begin(), dispDofs.end());
  Dofs.insert(Dofs.end(), piezoDofs.begin(), piezoDofs.end());

  indexType numDispDofs = static_cast<indexType>(dispDofs.size());
  indexType numPiezoDofs = static_cast<indexType>(piezoDofs.size());
  indexType numDofs = static_cast<indexType>(Dofs.size());
  

  elem->getSolution(pointers,Dofs,solution);

  Types::MatrixXX<prec> Bmat;
  Bmat.resize(9, numDofs);
  Bmat.setZero();
  std::vector<prec> xsi, eta, zeta, weight;
  elem->getGaussPoints(this->intOrderDisp + 1, weight, xsi, eta, zeta);
  for (auto i = 0; i < xsi.size(); ++i) {
    elem->getJacobian(*this->ptrCol, jacobi, xsi[i], eta[i], zeta[i]);
    elem->getH1Shapes(*this->ptrCol, this->intOrderDisp, jacobi, shapeDisp,
                      shapeDerivDisp, xsi[i], eta[i], zeta[i]);

    elem->getH1Shapes(*this->ptrCol, this->intOrderPiezo, jacobi, shapePiezo,
                      shapeDerivPiezo, xsi[i], eta[i], zeta[i]);

    for (auto l = 0; l < numDispDofs / 3; l++) {
      Bmat(0, l * 3) = shapeDerivDisp(0, l);
      Bmat(1, l * 3 + 1) = shapeDerivDisp(1, l);
      Bmat(2, l * 3 + 2) = shapeDerivDisp(2, l);

      Bmat(3, l * 3) = shapeDerivDisp(1, l);
      Bmat(3, l * 3 + 1) = shapeDerivDisp(0, l);

      Bmat(4, l * 3) = shapeDerivDisp(2, l);
      Bmat(4, l * 3 + 2) = shapeDerivDisp(0, l);

      Bmat(5, l * 3 + 1) = shapeDerivDisp(2, l);
      Bmat(5, l * 3 + 2) = shapeDerivDisp(1, l);
    }
    for (auto l = 0; l < numPiezoDofs / 3; l++) {
      Bmat(6, numDispDofs + l * 3) = shapeDerivPiezo(0, l);
      Bmat(7, numDispDofs + l * 3) = shapeDerivPiezo(1, l);
      Bmat(8, numDispDofs + l * 3) = shapeDerivPiezo(2, l);
    }
    epsilon = Bmat * solution;
    elem->getMaterialFormulation()->getMaterialData(epsilon, sigma, material);
    elem->projectOnVertsParaview(pointers,paraviewAdapter,gpDetJDummy,xsi[i], eta[i], zeta[i],weight[i],paraviewNames::weightName());
    elem->projectOnVertsParaview(pointers,paraviewAdapter,epsilon,xsi[i], eta[i], zeta[i],weight[i],paraviewNames::strainName());
    elem->projectOnVertsParaview(pointers,paraviewAdapter,sigma,xsi[i], eta[i], zeta[i],weight[i],paraviewNames::stressName());

  }
}

} // namespace HierAMuS

instantiate(Elementformulations::Hu_Washizu3DBeam)
