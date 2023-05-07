// Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
// SPDX-FileCopyrightText: 2023 Simon Klarmann
//
// SPDX-License-Identifier: BSD-3-Clause





#include "finiteElements/FaceConstraint.h"
#include "forwarddeclaration.h"
#include "geometry/Base.h"
#include "geometry/Edges.h"
#include "plot/vtkplotClassBase.h"
#include "shapefunctions/IntegrationsPoints/IntegrationPoints.h"
#include <iostream>


#include <elementFormulations/EL207_FaceConstraint.h>

#include <control/HandlingStructs.h>
#include <control/OutputHandler.h>
#include <pointercollection/pointercollection.h>

#include <equations/DegreeOfFreedom.h>
#include <equations/GenericNodes.h>

#include <finiteElements/GenericFiniteElement.h>
#include "finiteElements/Face.h"

#include <materials/GenericMaterialFormulation.h>

#include "math/geometry.h"


#include <elementFormulations/GenericElementFormulation.h>
#include <solver/GenericSolutionState.h>

#include <types/MatrixTypes.h>

#include <Eigen/Eigenvalues>
#include <vector>

#include <vtkCellType.h>

#include <Timer.h>

namespace HierAMuS::Elementformulations {

EL207_FaceConstraint::EL207_FaceConstraint(PointerCollection *ptrCol)
    : GenericElementFormulation(ptrCol) {}

EL207_FaceConstraint::~EL207_FaceConstraint() = default;

void EL207_FaceConstraint::readData(PointerCollection &pointers, ParameterList &list) {

  this->meshIdDisp = list.getIndexVal("meshIdDisp");
  this->meshIdRot = list.getIndexVal("meshIdRot");
  this->meshIdLam = list.getIndexVal("meshIdLam");
  this->meshIdMu = list.getIndexVal("meshIdMu");
  this->mode = list.getIndexVal("mode");
  this->intOrderDisp = list.getIndexVal("dispOrder");
  this->k = list.getPrecVal("k");

  auto Log = pointers.getSPDLogger();

  Log.info("\n{:-<100}\n"
                "*   Element 207, specified Options\n"
                "    Mesh id for displacement nodes:         {:>12}\n"
                "    Mesh id for rotation nodes:             {:>12}\n"
                "    Mesh id for disp lagrange parameters:   {:>12}\n"
                "    Mesh id for rot lagrange parameters:    {:>12}\n"
                "    Shape function order for displacements: {:>12}\n"
                "    Mode:                                   {:>12}\n"
                "{:-<100}\n",
                "",
                this->meshIdDisp,
                this->meshIdRot,
                this->meshIdLam,
                this->meshIdLam,
                this->intOrderDisp,
                this->mode,
                "");

  this->messageUnprocessed(pointers, list, "EL207_FaceConstraint");
}

void EL207_FaceConstraint::setDegreesOfFreedom(
  PointerCollection& pointers, FiniteElement::GenericFiniteElement *elem) {

  FiniteElement::FaceConstraint *facelem = static_cast<FiniteElement::FaceConstraint *>(elem);

  facelem->setH1Shapes(pointers, this->meshIdDisp, this->intOrderDisp);
  facelem->setVertexNodes(pointers, meshIdDisp);
  facelem->setVertexNodes(pointers, meshIdRot);
  facelem->setVertexNodes(pointers, meshIdLam);
  facelem->setVertexNodes(pointers, meshIdMu);
  
  
}

void EL207_FaceConstraint::AdditionalOperations(
  PointerCollection& pointers, FiniteElement::GenericFiniteElement *elem) {
  GenericElementFormulation::AdditionalOperations(pointers, elem);

  
}

auto EL207_FaceConstraint::getDofs(PointerCollection& pointers, FiniteElement::GenericFiniteElement *elem) -> std::vector<DegreeOfFreedom *>
{
  std::vector<DegreeOfFreedom*> Dofs;
  FiniteElement::FaceConstraint *facelem = static_cast<FiniteElement::FaceConstraint *>(elem);
  facelem->getH1Dofs(pointers, Dofs, this->meshIdDisp, this->intOrderDisp);
  facelem->getVertexDofs(pointers, Dofs, meshIdDisp);
  facelem->getVertexDofs(pointers, Dofs, meshIdRot);
  facelem->getVertexDofs(pointers, Dofs, meshIdLam);
  facelem->getVertexDofs(pointers, Dofs, meshIdMu);
  return Dofs;
}

void EL207_FaceConstraint::setTangentResidual(
  PointerCollection& pointers, FiniteElement::GenericFiniteElement *elem,
  Types::MatrixXX<prec> &stiffness, Types::VectorX<prec> &residual, std::vector<DegreeOfFreedom *> &Dofs) {

  switch (this->mode) {
  case 1: // Displacement based formulation
    this->setTangentResidualDispFormulation(pointers, elem, stiffness, residual,
                                            Dofs);
    break;
  default:
    throw std::runtime_error("Element formulation 207 called with wrong mode!");
  }
}


void EL207_FaceConstraint::toParaviewAdaper(PointerCollection &pointers,
                                     FiniteElement::GenericFiniteElement *elem,
                                     vtkPlotInterface &paraviewAdapter,
                                     ParaviewSwitch control) {

  int matNum = static_cast<int>(elem->getMaterial()->getNumber());
  switch (control) {
  case ParaviewSwitch::Mesh: {
    elem->geometryToParaview(pointers, paraviewAdapter, 0, matNum);

  } break;
  case ParaviewSwitch::Solution: {
    elem->H1SolutionToParaview(pointers, paraviewAdapter, 0, matNum,
                               this->meshIdDisp, this->intOrderDisp,
                               paraviewNames::DisplacementName());
  } break;
  case ParaviewSwitch::Weights: {
    elem->computeWeightsParaview(pointers, paraviewAdapter, 0, matNum);
  } break;
  case ParaviewSwitch::ProjectedValues: {
    
  }break;
  default:
  break;
  }
}

auto EL207_FaceConstraint::getHistoryDataStructure() -> const HistoryDataStructure & {
  return m_HistoryDataStructure;
}

auto EL207_FaceConstraint::getNumberOfIntergrationPoints(
  PointerCollection& pointers, FiniteElement::GenericFiniteElement *elem) -> indexType {
  auto GP = this->getIntegrationPoints(pointers, elem);
  return GP.getTotalGP();
}

void EL207_FaceConstraint::setTangentResidualDispFormulation(
  PointerCollection& pointers,
  FiniteElement::GenericFiniteElement *elem,
  Eigen::Matrix<prec, Eigen::Dynamic, Eigen::Dynamic> &stiffness,
  Eigen::Matrix<prec, Eigen::Dynamic, 1> &residual, std::vector<DegreeOfFreedom *> &Dofs) {


  FiniteElement::FaceConstraint *facelem = static_cast<FiniteElement::FaceConstraint *>(elem);


  Types::VectorX<prec> solution;
  Materials::MaterialTransferData materialData;

  Dofs.clear();
  facelem->getH1Dofs(pointers, Dofs, this->meshIdDisp, this->intOrderDisp);
  facelem->getVertexDofs(pointers, Dofs, meshIdDisp);
  facelem->getVertexDofs(pointers, Dofs, meshIdRot);
  facelem->getVertexDofs(pointers, Dofs, meshIdLam);
  facelem->getVertexDofs(pointers, Dofs, meshIdMu);

  indexType numDofs = static_cast<indexType>(Dofs.size());
  indexType numDispDofs = numDofs - 12;

  stiffness.resize(numDofs, numDofs);
  stiffness.setZero();

  residual.resize(numDofs);
  residual.setZero();

  elem->getSolution(pointers, Dofs, solution);

  auto GP = this->getIntegrationPoints(pointers, elem);
  GP.setOrder(intOrderDisp*intOrderDisp);

  auto Hist = elem->getHistoryDataIterator(pointers);
  Types::MatrixXX<prec> jaco;
  jaco.resize(2, 2);
  jaco = Types::MatrixXX<prec>::Identity(2, 2);
  Types::Vector3<prec> pointCoor = facelem->getVertexCoordinates(pointers);
  //prec sumDa = 0;
  Types::Matrix3X<prec> NMat;
  NMat.resize(3, numDispDofs);
  NMat.setZero();
  for (auto i : GP) {
    auto shapes = elem->getH1Shapes(pointers, intOrderDisp, jaco, i);
    for(auto j=0;j<shapes.shapes.rows();++j){
      NMat(0, 3*j+0) = shapes.shapes(j);
      NMat(1, 3*j+1) = shapes.shapes(j);
      NMat(2, 3*j+2) = shapes.shapes(j);
    }


    Types::Vector3<prec> surfaceCoor = facelem->getFaceCoordinates(pointers, i);
    Types::Vector3<prec> r0 = surfaceCoor - pointCoor;
    Types::Matrix33<prec> skewR = Math::Geometry::skewMatrix(r0);

    Types::Vector3<prec> G1 = facelem->getTangentG1(pointers, i);
    Types::Vector3<prec> G2 = facelem->getTangentG2(pointers, i);
    prec detJ = (G1.cross(G2)).norm();

    // dus*lam
    stiffness.block(0, numDispDofs+6,numDispDofs,3) += NMat.transpose() * detJ * i.weight;
    stiffness.block( numDispDofs+6, 0, 3, numDispDofs) += NMat * detJ * i.weight;

    // dus*skew(r) * mu
    stiffness.block(0, numDispDofs+9,numDispDofs,3) -= NMat.transpose() * skewR * detJ * i.weight;
    stiffness.block( numDispDofs+9, 0, 3, numDispDofs) -= skewR.transpose()*NMat * detJ * i.weight;

    // dub*lam
    stiffness.block(numDispDofs, numDispDofs+6, 3, 3) -= detJ * i.weight * Types::Matrix33<prec>::Identity(3,3);
    stiffness.block(numDispDofs+6, numDispDofs, 3, 3) -= detJ * i.weight * Types::Matrix33<prec>::Identity(3,3);

    // dub*skew(r) * mu
    stiffness.block(numDispDofs, numDispDofs+9, 3, 3) -= skewR * detJ * i.weight;
    stiffness.block(numDispDofs+9, numDispDofs, 3, 3) -= skewR.transpose() * detJ * i.weight;

    // dbeta*skew(r) * lam
    stiffness.block(numDispDofs+3, numDispDofs+6, 3, 3) -= skewR * detJ * i.weight;
    stiffness.block(numDispDofs+6, numDispDofs+3, 3, 3) += skewR * detJ * i.weight;

    // dbeta*skew(r).transpose()*skew(r) * mu
    stiffness.block(numDispDofs+3, numDispDofs+9, 3, 3) -= skewR.transpose()*skewR * detJ * i.weight;
    stiffness.block(numDispDofs+9, numDispDofs+3, 3, 3) -= skewR*skewR.transpose() * detJ * i.weight;

    //sumDa += detJ * i.weight;
    Hist.next();
    //sumDa += da;
  }
  stiffness*=this->k;
  residual = stiffness * solution;
  //std::cout << "sumDa: " << sumDa << std::endl;
  //std::cout << stiffness.eigenvalues() << std::endl;


}




auto EL207_FaceConstraint::getIntegrationPoints(PointerCollection& pointers, FiniteElement::GenericFiniteElement *elem)->IntegrationPoints{

  auto GP = elem->getIntegrationPoints(pointers);
  GP.setOrder(this->intOrderDisp*2);

  return GP;
}


const HistoryDataStructure EL207_FaceConstraint::m_HistoryDataStructure({{},{}},{{},{}});

} // namespace HierAMuS::Elementformulations
