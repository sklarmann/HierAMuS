// Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause



#include <elementFormulations/EL300_Solid3DLinear.h>
#include <pointercollection/pointercollection.h>

#include <control/HandlingStructs.h>
#include <control/OutputHandler.h>

#include <equations/DegreeOfFreedom.h>
#include <equations/GenericNodes.h>
#include <equations/NodeSet.h>
#include <equations/Nodetypes.h>

#include <finiteElements/GenericFiniteElement.h>

#include <geometry/Base.h>
#include <geometry/Volumes.h>


#include <materials/GenericMaterialFormulation.h>

#include "shapefunctions/IntegrationsPoints/dataClasses/GaussPoints.h"

#include <elementFormulations/GenericElementFormulation.h>
#include <solver/GenericSolutionState.h>

#include <types/MatrixTypes.h>

#include <Eigen/Eigenvalues>

namespace HierAMuS {
namespace Elementformulations {

EL300_Solid3DLinear::EL300_Solid3DLinear(PointerCollection *ptrCol)
    : GenericElementFormulation(ptrCol) {}

EL300_Solid3DLinear::~EL300_Solid3DLinear() {}

void EL300_Solid3DLinear::readData(PointerCollection &pointers,
                                   ParameterList &list) {
  this->meshIdDisp = list.getIndexVal("meshiddisp");
  this->intOrderDisp = list.getIndexVal("disporder");
  this->mode = list.getIndexVal("mode");
  
  auto Log = pointers.getSPDLogger();

  Log.info("\n{:-<100}\n"
                "*   Element 300, specified Options\n"
                "    Mesh id for displacement nodes:         {:>12}\n"
                "    Shape function order for displacements: {:>12}\n"
                "{:-<100}\n",
                "",
                this->meshIdDisp,
                this->intOrderDisp,
                "");


  this->messageUnprocessed(pointers, list, "EL300_Solid3DLinear");
}

void EL300_Solid3DLinear::setDegreesOfFreedom(
  PointerCollection& pointers, FiniteElement::GenericFiniteElement *elem) {

  elem->setH1Shapes(pointers, this->meshIdDisp, this->intOrderDisp);
}

void EL300_Solid3DLinear::AdditionalOperations(
  PointerCollection& pointers, FiniteElement::GenericFiniteElement *elem) {
  auto vol = elem->getVolume(pointers, 0);
  std::vector<indexType> faceNums;
  vol->getFaces(faceNums);
  for (auto i : faceNums) {
    auto face = pointers.getGeometryData()->getFace(i);
    auto nodes = face->getH1NodesInternal(pointers, this->meshIdDisp,
                                          this->intOrderDisp);
  }
}

auto EL300_Solid3DLinear::getDofs(PointerCollection& pointers, FiniteElement::GenericFiniteElement *elem) -> std::vector<DegreeOfFreedom*>{
  std::vector<DegreeOfFreedom *> ret;

  elem->getH1Dofs(pointers, ret, this->meshIdDisp, this->intOrderDisp);
  return ret;
}

void EL300_Solid3DLinear::setTangentResidual(
  PointerCollection& pointers, FiniteElement::GenericFiniteElement *elem,
  Types::MatrixXX<prec> &stiffness, Types::VectorX<prec> &residual, std::vector<DegreeOfFreedom *> &Dofs) {

  switch (this->mode) {
  case 1:
    this->setTangentResidualLinear(pointers, elem, stiffness, residual, Dofs);
    break;
  case 2:
    this->setTangentResidualNonLinear(pointers, elem, stiffness, residual,
                                      Dofs);
    break;
  }
}

void EL300_Solid3DLinear::setTangentResidualLinear(
  PointerCollection& pointers,
  FiniteElement::GenericFiniteElement *elem,
  Eigen::Matrix<prec, Eigen::Dynamic, Eigen::Dynamic> &stiffness,
  Eigen::Matrix<prec, Eigen::Dynamic, 1> &residual, std::vector<DegreeOfFreedom *> &Dofs) {

  Types::MatrixXX<prec> Bmat;
  Types::VectorX<prec> solution;

  Dofs.clear();
  elem->getH1Dofs(pointers, Dofs, this->meshIdDisp, this->intOrderDisp);

  indexType numDofs = static_cast<indexType>(Dofs.size());

  stiffness.resize(numDofs, numDofs);
  stiffness.setZero();

  residual.resize(numDofs);
  residual.setZero();

  elem->getSolution(pointers, Dofs, solution);

  Bmat.resize(6, numDofs);
  Bmat.setZero();

  auto GP = elem->getIntegrationPoints(pointers);
  GP.setOrder(this->intOrderDisp * 2);

  auto Hist = elem->getHistoryDataIterator(pointers);

  Materials::MaterialTransferData materialData;
  materialData.historyData = &Hist;
  //indexType cc = 0;
  for (auto i : GP) {
    auto jaco = elem->getJacobian(pointers, i);
    auto shapes = elem->getH1Shapes(pointers, this->intOrderDisp, jaco, i);

    auto Bmat = this->getLinearBMatrix(shapes);

    materialData.strains = Bmat * solution;
    elem->getMaterialFormulation(pointers)->getMaterialData(pointers, materialData,i);

    auto da = i.weight * jaco.determinant();
    stiffness += Bmat.transpose() * materialData.materialTangent * Bmat * da;
    residual += Bmat.transpose() * materialData.stresses * da;

    Hist.next();
  }
}

void EL300_Solid3DLinear::setTangentResidualNonLinear(
  PointerCollection& pointers,
  FiniteElement::GenericFiniteElement *elem,
  Eigen::Matrix<prec, Eigen::Dynamic, Eigen::Dynamic> &stiffness,
  Eigen::Matrix<prec, Eigen::Dynamic, 1> &residual, std::vector<DegreeOfFreedom *> &Dofs) {

  Dofs.clear();
  elem->getH1Dofs(pointers, Dofs, this->meshIdDisp, this->intOrderDisp);
  auto solution = elem->getSolution(pointers, Dofs);

  stiffness.resize(Dofs.size(), Dofs.size());
  residual.resize(Dofs.size());
  stiffness.setZero();
  residual.setZero();

  auto GP = elem->getIntegrationPoints(pointers);
  GP.setOrder(this->intOrderDisp * 2);

  auto Hist = elem->getHistoryDataIterator(pointers);

  Materials::MaterialTransferData materialData;
  materialData.historyData = &Hist;
  materialData.strains.resize(6);
  for (auto i : GP) {
    auto jaco = elem->getJacobian(pointers, i);
    auto shapes = elem->getH1Shapes(pointers, this->intOrderDisp, jaco, i);
    auto F = this->getDeformationGradient(shapes, solution);
    auto E =
        (F.transpose() * F - Types::Matrix33<prec>::Identity()) * prec(0.5);
    for (auto i = 0; i < 3; ++i) {
      materialData.strains(i) = E(i, i);
    }
    materialData.strains(3) = E(0, 1) * prec(2);
    materialData.strains(4) = E(0, 2) * prec(2);
    materialData.strains(5) = E(1, 2) * prec(2);
    elem->getMaterialFormulation(pointers)->getMaterialData(pointers, materialData,i);
    auto BMatrix = this->getNonLinearBMatrix(shapes, F);

    auto da = i.weight * jaco.determinant();
    stiffness +=
        BMatrix.transpose() * materialData.materialTangent * BMatrix * da;
    stiffness += this->getGeometricMatrix(shapes, materialData.stresses) * da;
    residual += BMatrix.transpose() * materialData.stresses * da;

    Hist.next();
  }
}

auto EL300_Solid3DLinear::getDeformationGradient(Geometry::H1Shapes &shapes,
                                                 Types::VectorX<prec> &solution)
    -> Types::Matrix33<prec> {
  Types::Matrix33<prec> F;
  F.setIdentity();

  for (auto i = 0; i < shapes.shapes.rows(); ++i) {
    for (auto j = 0; j < 3; ++j) {
      F(0, j) += solution(i * 3) * shapes.shapeDeriv(j, i);
      F(1, j) += solution(i * 3 + 1) * shapes.shapeDeriv(j, i);
      F(2, j) += solution(i * 3 + 2) * shapes.shapeDeriv(j, i);
    }
  }
  return F;
}

auto EL300_Solid3DLinear::getNonLinearBMatrix(Geometry::H1Shapes &shapes,
                                              Types::Matrix33<prec> &F)
    -> Types::Matrix6X<prec> {
  indexType nshapes = shapes.shapes.rows();
  Types::Matrix6X<prec> Bmat(6,nshapes*3);
  //Bmat.resize(6, nshapes * 3);
  Bmat.setZero();
  Types::Matrix33<prec> FT = F.transpose();
  auto dx1 = FT.block(0, 0, 1, 3);
  auto dx2 = FT.block(1, 0, 1, 3);
  auto dx3 = FT.block(2, 0, 1, 3);
  for (auto i = 0; i < nshapes; ++i) {
    Bmat.block(0, 3 * i, 1, 3) = shapes.shapeDeriv(0, i) * dx1;
    Bmat.block(1, 3 * i, 1, 3) = shapes.shapeDeriv(1, i) * dx2;
    Bmat.block(2, 3 * i, 1, 3) = shapes.shapeDeriv(2, i) * dx3;

    Bmat.block(3, 3 * i, 1, 3) =
        (shapes.shapeDeriv(0, i) * dx2 + shapes.shapeDeriv(1, i) * dx1);
    Bmat.block(4, 3 * i, 1, 3) =
        (shapes.shapeDeriv(0, i) * dx3 + shapes.shapeDeriv(2, i) * dx1);
    Bmat.block(5, 3 * i, 1, 3) =
        (shapes.shapeDeriv(1, i) * dx3 + shapes.shapeDeriv(2, i) * dx2);
  }

  return Bmat;
}

auto EL300_Solid3DLinear::getLinearBMatrix(Geometry::H1Shapes &shapes)
    -> Types::Matrix6X<prec> {
  Types::Matrix6X<prec> Bmat;
  indexType nshapes = shapes.shapes.rows();
  Bmat.resize(6, nshapes * 3);
  Bmat.setZero();
  for (auto l = 0; l < nshapes; ++l) {
    Bmat(0, l * 3) = shapes.shapeDeriv(0, l);
    Bmat(1, l * 3 + 1) = shapes.shapeDeriv(1, l);
    Bmat(2, l * 3 + 2) = shapes.shapeDeriv(2, l);

    Bmat(3, l * 3) = shapes.shapeDeriv(1, l);
    Bmat(3, l * 3 + 1) = shapes.shapeDeriv(0, l);

    Bmat(4, l * 3) = shapes.shapeDeriv(2, l);
    Bmat(4, l * 3 + 2) = shapes.shapeDeriv(0, l);

    Bmat(5, l * 3 + 1) = shapes.shapeDeriv(2, l);
    Bmat(5, l * 3 + 2) = shapes.shapeDeriv(1, l);
  }

  return Bmat;
}

auto EL300_Solid3DLinear::getGeometricMatrix(Geometry::H1Shapes &shapes,
                                            Types::VectorX<prec> &stress)
    -> Types::MatrixXX<prec> {
  Types::MatrixXX<prec> G;
  indexType nshapes = shapes.shapes.rows();
  G.resize(nshapes * 3, nshapes * 3);
  G.setZero();

  for (auto i = 0; i < nshapes; ++i) {
    for (auto j = 0; j < nshapes; ++j) {
      prec gik;
      gik = stress(0) * shapes.shapeDeriv(0, i) * shapes.shapeDeriv(0, j);
      gik += stress(1) * shapes.shapeDeriv(1, i) * shapes.shapeDeriv(1, j);
      gik += stress(2) * shapes.shapeDeriv(2, i) * shapes.shapeDeriv(2, j);
      gik += stress(3) * (shapes.shapeDeriv(0, i) * shapes.shapeDeriv(1, j) +
                          shapes.shapeDeriv(1, i) * shapes.shapeDeriv(0, j));
      gik += stress(4) * (shapes.shapeDeriv(0, i) * shapes.shapeDeriv(2, j) +
                          shapes.shapeDeriv(2, i) * shapes.shapeDeriv(0, j));
      gik += stress(5) * (shapes.shapeDeriv(1, i) * shapes.shapeDeriv(2, j) +
                          shapes.shapeDeriv(2, i) * shapes.shapeDeriv(1, j));

      G.block(i * 3, j * 3, 3, 3) = gik * Types::Matrix33<prec>::Identity();
    }
  }

  return G;
}

auto EL300_Solid3DLinear::getBTCBLinear(Geometry::H1Shapes &shapes,Types::MatrixXX<prec> &C, indexType nodeI, indexType nodeJ) -> Types::Matrix33<prec>
{

  Types::Matrix33<prec> res;

  prec N1x = shapes.shapeDeriv(0, nodeI);
  prec N1y = shapes.shapeDeriv(1, nodeI);
  prec N1z = shapes.shapeDeriv(2, nodeI);

  prec N2x = shapes.shapeDeriv(0, nodeJ);
  prec N2y = shapes.shapeDeriv(1, nodeJ);
  prec N2z = shapes.shapeDeriv(2, nodeJ);

  res(0,0) = N2x*(C(0,0)*N1x + C(0,3)*N1y + C(0,4)*N1z) + N2y*(C(0,3)*N1x + C(3,3)*N1y + C(3,4)*N1z) + N2z*(C(0,4)*N1x + C(3,4)*N1y + C(4,4)*N1z);
  res(0,1) = N2x*(C(0,3)*N1x + C(3,3)*N1y + C(3,4)*N1z) + N2y*(C(0,1)*N1x + C(1,3)*N1y + C(1,4)*N1z) + N2z*(C(0,5)*N1x + C(3,5)*N1y + C(4,5)*N1z);
  res(0,2) = N2x*(C(0,5)*N1x + C(3,5)*N1y + C(4,5)*N1z) + N2y*(C(0,4)*N1x + C(3,4)*N1y + C(4,4)*N1z) + N2z*(C(0,2)*N1x + C(2,3)*N1y + C(2,4)*N1z);
  res(1,0) = N2x*(C(0,1)*N1y + C(0,3)*N1x + C(0,5)*N1z) + N2y*(C(1,3)*N1y + C(3,3)*N1x + C(3,5)*N1z) + N2z*(C(1,4)*N1y + C(3,4)*N1x + C(4,5)*N1z);
  res(1,1) = N2x*(C(1,3)*N1y + C(3,3)*N1x + C(3,5)*N1z) + N2y*(C(1,1)*N1y + C(1,3)*N1x + C(1,5)*N1z) + N2z*(C(1,5)*N1y + C(3,5)*N1x + C(5,5)*N1z);
  res(1,2) = N2x*(C(1,5)*N1y + C(3,5)*N1x + C(5,5)*N1z) + N2y*(C(1,4)*N1y + C(3,4)*N1x + C(4,5)*N1z) + N2z*(C(1,2)*N1y + C(2,3)*N1x + C(2,5)*N1z);
  res(2,0) = N2x*(C(0,2)*N1z + C(0,4)*N1y + C(0,5)*N1x) + N2y*(C(2,3)*N1z + C(3,4)*N1y + C(3,5)*N1x) + N2z*(C(2,4)*N1z + C(4,4)*N1y + C(4,5)*N1x);
  res(2,1) = N2x*(C(2,3)*N1z + C(3,4)*N1y + C(3,5)*N1x) + N2y*(C(1,2)*N1z + C(1,4)*N1y + C(1,5)*N1x) + N2z*(C(2,5)*N1z + C(4,5)*N1y + C(5,5)*N1x);
  res(2,2) = N2x*(C(2,5)*N1z + C(4,5)*N1y + C(5,5)*N1x) + N2y*(C(2,4)*N1z + C(4,4)*N1y + C(4,5)*N1x) + N2z*(C(2,2)*N1z + C(2,4)*N1y + C(2,5)*N1x);

  return res;
}

void EL300_Solid3DLinear::toParaviewAdaper(
    PointerCollection &pointers, FiniteElement::GenericFiniteElement *elem,
    vtkPlotInterface &paraviewAdapter, ParaviewSwitch control) {
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
    auto GP = elem->getIntegrationPoints(pointers);
    GP.setOrder(this->intOrderDisp * 2);

    auto Hist = elem->getHistoryDataIterator(pointers);

    Materials::MaterialTransferData materialData;
    materialData.historyData = &Hist;
    materialData.strains.resize(6);
    std::vector<DegreeOfFreedom *> Dofs;

    elem->getH1Dofs(pointers, Dofs, this->meshIdDisp, this->intOrderDisp);
    auto solution = elem->getSolution(pointers, Dofs);

    for (auto i : GP) {
      auto jaco = elem->getJacobian(pointers, i);
      auto shapes = elem->getH1Shapes(pointers, this->intOrderDisp, jaco, i);

      if (this->mode == 1) {
        auto BMat = this->getLinearBMatrix(shapes);
        materialData.strains = BMat * solution;
      } else {
        auto F = this->getDeformationGradient(shapes, solution);
        auto E =
            (F.transpose() * F - Types::Matrix33<prec>::Identity()) * prec(0.5);
        materialData.strains(0) = E(0, 0);
        materialData.strains(1) = E(1, 1);
        materialData.strains(2) = E(2, 2);
        materialData.strains(3) = E(0, 1);
        materialData.strains(4) = E(0, 2);
        materialData.strains(5) = E(1, 2);
      }

      elem->getMaterialFormulation(pointers)->getMaterialData(pointers, materialData,i);
      elem->projectDataToParaviewVertices(
          pointers, paraviewAdapter, 0, matNum, this->intOrderDisp, i,
          materialData.strains, 6, paraviewNames::strainName());
      elem->projectDataToParaviewVertices(
          pointers, paraviewAdapter, 0, matNum, this->intOrderDisp, i,
          materialData.stresses, 6, paraviewNames::stressName());
      Hist.next();
    }
  }break;
  default:
    break;
  }
}

auto EL300_Solid3DLinear::getHistoryDataStructure()
    -> const HistoryDataStructure & {
  return m_HistoryDataStructure;
}

auto EL300_Solid3DLinear::getNumberOfIntergrationPoints(
  PointerCollection& pointers, FiniteElement::GenericFiniteElement *elem) -> indexType {
  auto GP = elem->getIntegrationPoints(pointers);
  GP.setOrder(this->intOrderDisp * 2);
  return GP.getTotalGP();
}

const HistoryDataStructure EL300_Solid3DLinear::m_HistoryDataStructure({},{});

} // namespace Elementformulations
} // namespace HierAMuS
