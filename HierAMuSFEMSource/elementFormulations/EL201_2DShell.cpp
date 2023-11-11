// Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause

#include "elementFormulations/EL201_2DShell.h"

#include "geometry/GeometryBaseData.h"
#include "geometry/Edges/EdgesData.h"


#include "finiteElements/Face.h"
#include <finiteElements/GenericFiniteElement.h>

#include <materials/GenericMaterialFormulation.h>

#include "control/ParameterList.h"

#include <elementFormulations/GenericElementFormulation.h>
#include <solver/GenericSolutionState.h>

#include <types/MatrixTypes.h>

#include <Eigen/Eigenvalues>
#include <vector>

#include <vtkCellType.h>

#include <Timer.h>

namespace HierAMuS::Elementformulations {

EL201_2DShell::EL201_2DShell(PointerCollection *ptrCol)
    : GenericElementFormulationInterface(ptrCol) {}

EL201_2DShell::~EL201_2DShell() = default;

void EL201_2DShell::readData(PointerCollection &pointers, ParameterList &list) {

  this->meshIdDisp = list.getIndexVal("meshiddisp");
  this->mode = list.getIndexVal("mode");
  this->intOrderDisp = list.getIndexVal("disporder");
  this->meshIdStre = 2;
  this->intOrderStre = intOrderDisp - 1;

  auto &Logger = pointers.getSPDLogger();

  Logger.info("\n{:-<100}\n"
              "*   Element 201, specified Options\n"
              "    Mesh id for displacement nodes:          {:>12}\n"
              "    Shape function order for displacements:  {:>12}\n"
              "\n{:-<100}\n",
              "", this->meshIdDisp, this->intOrderDisp, "");

  this->messageUnprocessed(pointers, list, "EL201_2DShell");
}

void EL201_2DShell::setDegreesOfFreedom(PointerCollection &pointers,
                                        FiniteElement::Face &elem) {

  elem.setH1Shapes(pointers, this->meshIdDisp, this->intOrderDisp);
  auto GP = elem.getIntegrationPoints(pointers);
  GP.setOrder(this->intOrderDisp * 2);

  if (this->mode == 4) {
    elem.setL2Shapes(pointers, this->meshIdStre, this->intOrderStre);
  }
}

void EL201_2DShell::AdditionalOperations(PointerCollection &pointers,
                                         FiniteElement::Face &elem) {

  auto ip = this->getIntegrationPoints(pointers, &elem);
  for (auto gp : ip) {
    elem.getMaterialFormulation(pointers)->initRVE(pointers, gp);
  }

  elem.setAllNodeBoundaryConditionMeshId(pointers, this->meshIdDisp, 2);
  auto GP = this->getIntegrationPoints(pointers, &elem);

  if (this->mode == 4) {
    elem.setAllNodeBoundaryConditionMeshId(pointers, this->meshIdStre, 1);
    elem.setAllNodeBoundaryConditionMeshId(pointers, this->meshIdStre, 2);
  }

  // if (this->mode == 2) {
  //   Types::MatrixXX<prec> jaco0;
  //   IntegrationPoint ip;
  //   ip.xi = prec(0);
  //   ip.eta = prec(0);
  //   jaco0 = elem->getJacobian(pointers, ip);
  //   auto hist = elem->getHistoryDataIterator(pointers);
  //   auto jacoHist = hist.getFieldElementConst(1);
  //   jacoHist = jaco0;
  //   auto centGrav = hist.getFieldElementConst(0);

  //  prec dA0 = prec(4) * jaco0.determinant();
  //  auto GP = this->getIntegrationPoints(pointers, elem);

  //  prec xiQ = prec(0), etaQ = prec(0);
  //  for (auto ip : GP) {
  //    Types::MatrixXX<prec> jac = elem->getJacobian(pointers, ip);
  //    prec dA = jac.determinant() * ip.weight;
  //    xiQ += ip.xi * dA;
  //    etaQ += ip.eta * dA;
  //  }
  //  centGrav(0) = xiQ / dA0;
  //  centGrav(1) = etaQ / dA0;
  //}
}

auto EL201_2DShell::getDofs(PointerCollection &pointers,
                            FiniteElement::GenericFiniteElement *elem)
    -> std::vector<DegreeOfFreedom *> {
  std::vector<DegreeOfFreedom *> Dofs;
  elem->getH1Dofs(pointers, Dofs, this->meshIdDisp, this->intOrderDisp);
  return Dofs;
}

void EL201_2DShell::setTangentResidual(PointerCollection &pointers,
                                       FiniteElement::Face &elem,
                                       Types::MatrixXX<prec> &stiffness,
                                       Types::VectorX<prec> &residual,
                                       std::vector<DegreeOfFreedom *> &Dofs) {

  switch (this->mode) {
  case 1: // Displacement based formulation
    this->setTangentResidualDispFormulation(pointers, elem, stiffness, residual,
                                            Dofs);
    break;
  case 2: // Hu-Washizu formulation
    this->setTangentResidualHuWashizuFormulation(pointers, elem, stiffness,
                                                 residual, Dofs);
    break;
  case 3:
    this->setTangentResidualDispNonLinear(pointers, elem, stiffness, residual,
                                          Dofs);
    break;
  case 4:
    this->setTangentResidualDispAndStreFormulation(pointers, elem, stiffness,
                                                   residual, Dofs);
    break;
  default:
    throw std::runtime_error("Element formulation 201 called with wrong mode!");
  }
}

void EL201_2DShell::toParaviewAdaper(PointerCollection &pointers,
                                     FiniteElement::Face &elem,
                                     vtkPlotInterface &paraviewAdapter,
                                     ParaviewSwitch control) {

  int matNum = static_cast<int>(elem.getMaterial()->getNumber());
  switch (control) {
  case ParaviewSwitch::Mesh: {
    elem.geometryToParaview(pointers, paraviewAdapter, 0, matNum);

  } break;
  case ParaviewSwitch::Solution: {
    elem.H1SolutionToParaview(pointers, paraviewAdapter, 0, matNum,
                              this->meshIdDisp, this->intOrderDisp,
                              paraviewNames::DisplacementName());
  } break;
  case ParaviewSwitch::Weights: {
    elem.computeWeightsParaview(pointers, paraviewAdapter, 0, matNum);
  } break;
  case ParaviewSwitch::ProjectedValues: {
    indexType matNum = elem.getMaterial()->getNumber();
    auto GP = this->getIntegrationPoints(pointers, &elem);
    auto Hist = elem.getHistoryDataIterator(pointers);

    Types::VectorX<prec> solution;
    Materials::MaterialTransferData materialData;
    std::vector<DegreeOfFreedom *> Dofs;

    elem.getH1Dofs(pointers, Dofs, this->meshIdDisp, this->intOrderDisp);
    elem.getSolution(pointers, Dofs, solution);
    Types::Matrix3X<prec> Bmat = Types::Matrix3X<prec>::Zero(3,Dofs.size());

    Types::VectorX<prec> eps = Types::VectorX<prec>::Zero(6);
    Types::VectorX<prec> sig = Types::VectorX<prec>::Zero(6);

    materialData.historyData = &Hist;
    for (auto i : GP) {
      auto jacobi = elem.getJacobian(pointers, i);
      auto shapes = elem.getH1Shapes(pointers, this->intOrderDisp, jacobi, i);
      Bmat = this->getBMatrix(shapes.shapeDeriv, Dofs.size());
      materialData.strains = Bmat * solution;
      auto internalMaterialData =
          elem.getMaterialFormulation(pointers)->getInternalVariables(
              pointers, materialData);
      elem.getMaterialFormulation(pointers)->getMaterialData(pointers,
                                                             materialData, i);
      eps.block(0, 0, 3, 1) = materialData.strains;
      sig.block(0, 0, 3, 1) = materialData.stresses;

      eps(3) = eps(2);
      sig(3) = sig(2);
      eps(2) = prec(0);
      sig(2) = prec(0);

      for (auto j : internalMaterialData) {
        elem.projectDataToParaviewVertices(pointers, paraviewAdapter, 0, matNum,
                                           this->intOrderDisp, i, j.second,
                                           j.second.rows(), j.first);
      }

      elem.projectDataToParaviewVertices(pointers, paraviewAdapter, 0, matNum,
                                         this->intOrderDisp, i, eps, 6,
                                         paraviewNames::strainName());
      elem.projectDataToParaviewVertices(pointers, paraviewAdapter, 0, matNum,
                                         this->intOrderDisp, i, sig, 6,
                                         paraviewNames::stressName());
      Hist.next();
    }
  } break;
  default:
    break;
  }
}

auto EL201_2DShell::getHistoryDataStructure() -> const HistoryDataStructure & {
  return m_HistoryDataStructure;
}

auto EL201_2DShell::getNumberOfIntergrationPoints(
    PointerCollection &pointers, FiniteElement::GenericFiniteElement *elem)
    -> indexType {
  auto GP = this->getIntegrationPoints(pointers, elem);
  return GP.getTotalGP();
}

void EL201_2DShell::setTangentResidualDispFormulation(
    PointerCollection &pointers, FiniteElement::Face &elem,
    Eigen::Matrix<prec, Eigen::Dynamic, Eigen::Dynamic> &stiffness,
    Eigen::Matrix<prec, Eigen::Dynamic, 1> &residual,
    std::vector<DegreeOfFreedom *> &Dofs) {

  Types::VectorX<prec> solution;
  Materials::MaterialTransferData materialData;

  Dofs.clear();
  elem.getH1Dofs(pointers, Dofs, this->meshIdDisp, this->intOrderDisp);

  auto numDofs = static_cast<indexType>(Dofs.size());

  stiffness.resize(numDofs, numDofs);
  stiffness.setZero();

  residual.resize(numDofs);
  residual.setZero();

  elem.getSolution(pointers, Dofs, solution);

  auto GP = this->getIntegrationPoints(pointers, &elem);

  auto Hist = elem.getHistoryDataIterator(pointers);

  // prec sumDa = 0;
  for (auto i : GP) {
    Types::Matrix22<prec> jacob = elem.getJacobian(pointers, i);
    auto shapes = elem.getH1Shapes(pointers, this->intOrderDisp, jacob, i);
    // elem->getH1Shapes(*this->ptrCol, this->intOrderDisp, jacob,
    // shapes.shapes, shapes.shapeDeriv, i);
    Types::Matrix3X<prec> Bmat = getBMatrix(shapes.shapeDeriv, numDofs);
    materialData.strains = Bmat * solution;
    materialData.historyData = &Hist;
    elem.getMaterialFormulation(pointers)->getMaterialData(pointers,
                                                           materialData, i);
    prec da = jacob.determinant() * i.weight;
    stiffness += Bmat.transpose() * materialData.materialTangent * Bmat * da;
    residual += Bmat.transpose() * materialData.stresses * da;

    Hist.next();
    // sumDa += da;
  }
  // std::cout << "sumDa: " << sumDa << std::endl;
  // std::cout << stiffness.eigenvalues() << std::endl;
}

Types::Matrix3X<prec>
EL201_2DShell::getBMatrix(const Types::Matrix2X<prec> &shapeDeriv,
                          indexType numDofs) {

  Types::Matrix3X<prec> Bmat = Types::Matrix3X<prec>::Zero(3,numDofs);
  for (auto j = 0; j < numDofs / 3; ++j) {
    Bmat(0, j * 3) = shapeDeriv(0, j);
    Bmat(1, j * 3 + 1) = shapeDeriv(1, j);
    Bmat(2, j * 3) = shapeDeriv(1, j);
    Bmat(2, j * 3 + 1) = shapeDeriv(0, j);
  }
  return Bmat;
}

void EL201_2DShell::setTangentResidualHuWashizuFormulation(
    PointerCollection &pointers, FiniteElement::Face &elem,
    Eigen::Matrix<prec, Eigen::Dynamic, Eigen::Dynamic> &stiffness,
    Eigen::Matrix<prec, Eigen::Dynamic, 1> &residual,
    std::vector<DegreeOfFreedom *> &Dofs) {

  Types::Matrix2X<prec> shapeDeriv;
  Types::VectorXT<prec> shape;
  Types::VectorX<prec> shapeT;
  Types::Matrix2X<prec> coorDeriv;
  Types::VectorXT<prec> coor;
  Types::VectorX<prec> solution;
  std::vector<GenericNodes *> dispNodes;
  Materials::MaterialTransferData materialData;

  Dofs.clear();
  elem.getH1Dofs(pointers, Dofs, this->meshIdDisp, this->intOrderDisp);

  auto numDofs = static_cast<indexType>(Dofs.size());

  stiffness.resize(numDofs, numDofs);
  stiffness.setZero();

  residual.resize(numDofs);
  residual.setZero();

  elem.getSolution(pointers, Dofs, solution);

  prec da;
  da = 0;

  Types::Matrix3X<prec> Bmat = Types::Matrix3X<prec>::Zero(3,numDofs);

  Types::Matrix3X<prec> Amat;
  Amat =
      this->getLocalStrainStressInterpolation(pointers, elem, prec(0), prec(0));
  indexType AmatParam = Amat.cols();
  // Amat = this->getLocalStrainStressInterpolation();

  // indexType approxGP = (this->intOrderDisp + 1) * (this->intOrderDisp + 1);

  // std::vector<prec> xsi, eta, weight;
  // xsi.reserve(approxGP), eta.reserve(approxGP), weight.reserve(approxGP);

  auto GP = this->getIntegrationPoints(pointers, &elem);

  Types::MatrixXX<prec> G, H, L;
  G.resize(numDofs, AmatParam);
  H.resize(AmatParam, AmatParam);
  L.resize(AmatParam, AmatParam);

  G.setZero();
  H.setZero();
  L.setZero();

  IntegrationPoint ip;
  ip.xi = prec(0);
  ip.eta = prec(0);
  auto jaco0 = elem.getJacobian(pointers, ip);

  auto &tEdge = elem.getEdge(pointers, 0);

  Types::Matrix22<prec> epsShape;

  for (auto &ip : GP) {
    epsShape.setZero();

    IntegrationPoint ippTemp;
    ippTemp.xi = ip.xi;
    auto ShapesXi = tEdge.getH1Shapes(indexType(1), ippTemp);
    ippTemp.xi = ip.eta;
    auto ShapesEta = tEdge.getH1Shapes(indexType(1), ippTemp);

    for (auto j = 0; j < 2; ++j) {
      epsShape(0, j) +=
          jaco0(0, 0) * ShapesXi.shapes(j) + jaco0(1, 0) * ShapesEta.shapes(j);
      epsShape(1, j) +=
          jaco0(0, 1) * ShapesXi.shapes(j) + jaco0(1, 1) * ShapesEta.shapes(j);
    }

    Amat =
        this->getLocalStrainStressInterpolation(pointers, elem, ip.xi, ip.eta);

    auto jacobi = elem.getJacobian(pointers, ip);
    auto elemH1Shapes =
        elem.getH1Shapes(pointers, this->intOrderDisp, jacobi, ip);

    Bmat = getBMatrix(shapeDeriv, numDofs);

    materialData.strains = Bmat * solution;

    elem.getMaterialFormulation(pointers)->getMaterialData(
        pointers, materialData, GP.getIntegrationPoint());

    da = (jacobi.determinant()) * ip.weight;
    G += Bmat.transpose() * Amat * da;
    H -= Amat.transpose() * Amat * da;
    L += Amat.transpose() * materialData.materialTangent * Amat * da;
  }

  stiffness = G * (H * L.inverse() * H.transpose()).inverse() * G.transpose();

  // Eigen::EigenSolver<Types::MatrixXX<prec>> es(stiffness);
  // std::cout << "Eigenvalues of J are:\n " << es.eigenvalues() << std::endl;

  residual = stiffness * solution;
}

Types::Matrix3X<prec> EL201_2DShell::getLocalStrainStressInterpolation(
    PointerCollection &pointers, FiniteElement::Face &elem, prec xi, prec eta) {

  auto Hist = elem.getHistoryDataIterator(pointers);

  auto centGrav = Hist.getFieldElementConst(0);
  auto jac0 = Hist.getFieldElementConst(1);

  Types::Matrix3X<prec> Amat;
  Amat.resize(3, 5);
  Amat.setZero();

  Amat(0, 0) = prec(1);
  Amat(0, 3) = jac0(0, 0) * jac0(0, 0) * (eta - centGrav(1));
  Amat(0, 4) = jac0(1, 0) * jac0(1, 0) * (xi - centGrav(0));

  Amat(1, 1) = prec(1);
  Amat(1, 3) = jac0(0, 1) * jac0(0, 1) * (eta - centGrav(1));
  Amat(1, 4) = jac0(1, 1) * jac0(1, 1) * (xi - centGrav(0));

  Amat(2, 2) = prec(1);
  Amat(2, 3) = jac0(0, 0) * jac0(0, 1) * (eta - centGrav(1));
  Amat(2, 4) = jac0(1, 0) * jac0(1, 1) * (xi - centGrav(0));

  return Amat;
}

void EL201_2DShell::setTangentResidualDispNonLinear(
    PointerCollection &pointers, FiniteElement::Face &elem,
    Eigen::Matrix<prec, Eigen::Dynamic, Eigen::Dynamic> &stiffness,
    Eigen::Matrix<prec, Eigen::Dynamic, 1> &residual,
    std::vector<DegreeOfFreedom *> &Dofs) {

  elem.getH1Dofs(pointers, Dofs, this->meshIdDisp, this->intOrderDisp);
  stiffness.resize(Dofs.size(), Dofs.size());
  residual.resize(Dofs.size());
  stiffness.setZero();
  residual.setZero();

  auto GP = this->getIntegrationPoints(pointers, &elem);
  auto Hist = elem.getHistoryDataIterator(pointers);

  Types::VectorX<prec> sol = elem.getSolution(pointers, Dofs);
  Materials::MaterialTransferData materialData;
  materialData.strains.resize(3);
  materialData.historyData = &Hist;
  for (auto i : GP) {
    Types::Matrix22<prec> jacobi = elem.getJacobian(pointers, i);
    auto shapes = elem.getH1Shapes(pointers, this->intOrderDisp, jacobi, i);
    Types::Matrix22<prec> F = this->getDeformationGradient(shapes, sol);
    Types::Matrix3X<prec> BMat = this->getNonLinearBMatrix(shapes, F);
    Types::Matrix22<prec> E =
        prec(0.5) * (F.transpose() * F - Types::Matrix22<prec>::Identity());
    materialData.strains(0) = E(0, 0);
    materialData.strains(1) = E(1, 1);
    materialData.strains(2) = prec(2) * E(1, 0);
    elem.getMaterialFormulation(pointers)->getMaterialData(pointers,
                                                           materialData, i);
    prec da = i.weight * jacobi.determinant();
    stiffness += BMat.transpose() * materialData.materialTangent * BMat * da;
    stiffness += this->getGeometricMatrix(shapes, materialData.stresses) * da;
    residual += BMat.transpose() * materialData.stresses * da;

    Hist.next();
  }
}
auto EL201_2DShell::getDeformationGradient(Geometry::H1Shapes &shapes,
                                           Types::VectorX<prec> &solution)
    -> Types::Matrix22<prec> {
  Types::Matrix22<prec> F;
  F = Types::Matrix22<prec>::Identity();
  for (auto i = 0; i < 2; ++i) {
    for (auto j = 0; j < 2; ++j) {
      for (auto k = 0; k < shapes.shapes.rows(); ++k) {
        F(i, j) += shapes.shapeDeriv(j, k) * solution(3 * k + i);
      }
    }
  }
  return F;
}

auto EL201_2DShell::getNonLinearBMatrix(Geometry::H1Shapes &shapes,
                                        Types::Matrix22<prec> &F)
    -> Types::Matrix3X<prec> {

  indexType nshapes = shapes.shapes.rows();
  Types::Matrix3X<prec> Bmat;
  Bmat.resize(3, nshapes * 3);
  Bmat.setZero();
  Types::Matrix22<prec> FT = F.transpose();

  for (auto i = 0; i < nshapes; ++i) {
    Bmat.block(0, 3 * i, 1, 2) = FT.block(0, 0, 1, 2) * shapes.shapeDeriv(0, i);
    Bmat.block(1, 3 * i, 1, 2) = FT.block(1, 0, 1, 2) * shapes.shapeDeriv(1, i);
    Bmat.block(2, 3 * i, 1, 2) =
        FT.block(1, 0, 1, 2) * shapes.shapeDeriv(0, i) +
        FT.block(0, 0, 1, 2) * shapes.shapeDeriv(1, i);
  }

  return Bmat;
}

auto EL201_2DShell::getGeometricMatrix(Geometry::H1Shapes &shapes,
                                       Types::VectorX<prec> &stress)
    -> Types::MatrixXX<prec> {
  indexType nshapes = shapes.shapes.rows();
  Types::MatrixXX<prec> G;
  G.resize(nshapes * 3, nshapes * 3);
  G.setZero();

  for (auto i = 0; i < nshapes; ++i) {
    for (auto j = 0; j < nshapes; ++j) {
      prec gik;
      gik = stress(0) * shapes.shapeDeriv(0, i) * shapes.shapeDeriv(0, j);
      gik += stress(1) * shapes.shapeDeriv(1, i) * shapes.shapeDeriv(1, j);
      gik += stress(2) * (shapes.shapeDeriv(0, i) * shapes.shapeDeriv(1, j) +
                          shapes.shapeDeriv(1, i) * shapes.shapeDeriv(0, j));

      G.block(i * 3, j * 3, 2, 2) = gik * Types::Matrix22<prec>::Identity();
    }
  }

  return G;
}

auto EL201_2DShell::getIntegrationPoints(
    PointerCollection &pointers, FiniteElement::GenericFiniteElement *elem)
    -> IntegrationPoints {

  auto GP = elem->getIntegrationPoints(pointers);
  GP.setOrder(this->intOrderDisp * 2);

  return GP;
}

void EL201_2DShell::setTangentResidualDispAndStreFormulation(
    PointerCollection &pointers, FiniteElement::Face &elem,
    Eigen::Matrix<prec, Eigen::Dynamic, Eigen::Dynamic> &stiffness,
    Eigen::Matrix<prec, Eigen::Dynamic, 1> &residual,
    std::vector<DegreeOfFreedom *> &Dofs) {

  Types::VectorX<prec> solution;
  Materials::MaterialTransferData materialData;
  Eigen::Matrix<prec, Eigen::Dynamic, Eigen::Dynamic> Kuu;
  Eigen::Matrix<prec, Eigen::Dynamic, Eigen::Dynamic> Kup;
  // Eigen::Matrix<prec, Eigen::Dynamic, Eigen::Dynamic> Kpu;
  Eigen::Matrix<prec, Eigen::Dynamic, Eigen::Dynamic> Kpp;
  indexType imagOrder;

  // Materialeigenschaften setzen
  prec E = 100;
  prec nu = 0.25;

  // umrechnen
  prec myu = E / (2 * (1 + nu));
  prec kappa = E / (3 * (1 - 2 * nu));
  Types::Matrix33<prec> Cdev;

  prec a = 4.0 / 3.0 * myu;
  prec b = -2.0 / 3.0 * myu;
  Cdev(0, 0) = a;
  Cdev(1, 0) = b;
  Cdev(2, 0) = 0;
  Cdev(0, 1) = b;
  Cdev(0, 2) = 0;
  Cdev(1, 1) = a;
  Cdev(2, 2) = myu;
  Cdev(2, 1) = 0;
  Cdev(1, 2) = 0;

  Dofs.clear();
  elem.getH1Dofs(pointers, Dofs, this->meshIdDisp, this->intOrderDisp);

  auto numDofsDisp = static_cast<indexType>(Dofs.size());

  elem.getL2Dofs(pointers, Dofs, this->meshIdStre, this->intOrderStre);
  auto numDofsStre = static_cast<indexType>(Dofs.size()) - numDofsDisp;

  Kuu.resize(numDofsDisp, numDofsDisp);
  Kuu.setZero();
  Kup.resize(numDofsDisp, numDofsStre);
  Kup.setZero();
  Kpp.resize(numDofsStre, numDofsStre);
  Kpp.setZero();

  residual.resize(numDofsStre + numDofsDisp);
  residual.setZero();

  elem.getSolution(pointers, Dofs, solution);

  auto GP = this->getIntegrationPoints(pointers, &elem);

  for (auto i : GP) {
    Types::Matrix22<prec> jacob = elem.getJacobian(pointers, i);
    auto H1shapes = elem.getH1Shapes(pointers, this->intOrderDisp, jacob, i);
    if (this->intOrderStre == 0) {
      imagOrder = 1;
    } else {
      imagOrder = this->intOrderStre;
    }
    auto L2shapes = elem.getL2Shapes(pointers, imagOrder, jacob, i);

    Types::Matrix3X<prec> Bmat = getBMatrix(H1shapes.shapeDeriv, numDofsDisp);
    Types::VectorXT<prec> Hmat;
    if (intOrderStre == 0) {
      Hmat.resize(1, 3);
      Hmat.setOnes();
    } else {
      Hmat = getHMatrix(L2shapes.shapes, numDofsStre);
    }

    prec da = jacob.determinant() * i.weight;
    Types::VectorXT<prec> Bvolmat;
    Bvolmat = getBvolMatrix(H1shapes.shapeDeriv, numDofsDisp);

    Kuu += Bmat.transpose() * Cdev * Bmat * da;
    Kup += Bvolmat.transpose() * Hmat * da;
    Kpp -= Hmat.transpose() * 1 / kappa * Hmat * da;
  }
  indexType row = numDofsDisp + numDofsStre;
  stiffness.resize(row, row);
  stiffness.block(0, 0, numDofsDisp, numDofsDisp) = Kuu;
  stiffness.block(0, numDofsDisp, numDofsDisp, numDofsStre) = Kup;
  stiffness.block(numDofsDisp, 0, numDofsStre, numDofsDisp) = Kup.transpose();
  stiffness.block(numDofsDisp, numDofsDisp, numDofsStre, numDofsStre) = Kpp;
  /*
  for (int i = 0; i < stiffness.rows(); ++i) {
      std::cout << "[";
      for (int j = 0; j < stiffness.cols(); ++j) {
          std::cout << stiffness(i, j) << ",\t";
          if ((j + 1) % 3 == 0) {
              std::cout << "\n";
          }
      }
      std::cout << "]\n";
  }
  */
  residual = stiffness * solution;
}
Types::VectorXT<prec>
EL201_2DShell::getHMatrix(const Types::VectorX<prec> &shapes,
                          indexType numDofs) {

  Types::VectorXT<prec> Hmat;
  Hmat.resize(1, numDofs);
  Hmat.setZero();
  for (auto j = 0; j < numDofs / 3; ++j) {
    Hmat(0, j * 3) = shapes(j);
    // Hmat(1, j * 3 + 1) = shapes(j);
  }
  return Hmat;
}

Types::VectorXT<prec>
EL201_2DShell::getBvolMatrix(const Types::Matrix2X<prec> &shapeDeriv,
                             indexType numDofs) {

  Types::VectorXT<prec> Bvolmat;
  Bvolmat.resize(1, numDofs);
  Bvolmat.setZero();
  for (auto j = 0; j < numDofs / 3; ++j) {
    Bvolmat(0, j * 3) = shapeDeriv(0, j);
    Bvolmat(0, j * 3 + 1) = shapeDeriv(1, j);
  }
  return Bvolmat;
}

const HistoryDataStructure
    EL201_2DShell::m_HistoryDataStructure({{2, 1}, {2, 2}}, {{1, 1}, {1, 1}});

} // namespace HierAMuS::Elementformulations
