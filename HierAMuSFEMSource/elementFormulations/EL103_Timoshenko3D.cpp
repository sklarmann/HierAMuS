// Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause

#include "plot/vtkplotClassBase.h"


#include "control/ParameterList.h"

#include <elementFormulations/EL103_Timoshenko3D.h>
#include <elementFormulations/GenericElementFormulation.h>

#include "shapefunctions/IntegrationsPoints/IntegrationPoints.h"

#include "finiteElements/Edge.h"
#include "materials/GenericMaterialFormulation.h"

#include <geometry/GeometryBaseData.h>
#include <geometry/VertexData.h>


#include "pointercollection/pointercollection.h"
#include <solver/GenericSolutionState.h>

#include <Eigen/Dense>

#include "PropfunctionHandler.h"

#include <stdexcept>
#include <types/MatrixTypes.h>
#include <vector>

#include <vtkCellType.h>

#include "spdlog/fmt/ostr.h"

namespace HierAMuS::Elementformulations {

EL103_Timoshenko3D::EL103_Timoshenko3D(PointerCollection *ptrCol)
    : GenericElementFormulationInterface(ptrCol) {}

EL103_Timoshenko3D::~EL103_Timoshenko3D() = default;

void EL103_Timoshenko3D::readData(PointerCollection &pointers,
                                  ParameterList &list) {

  m_RVECmatInit = false;
  this->geo = list.getIndexVal("geo");
  this->thermal = list.getIndexVal("thermal");
  this->fe2 = list.getIndexVal("fe2");

  if (thermal == 0) {
    if (this->fe2 == 0) {
      if (this->geo == 0) {
        this->mode = 1;
      } else if (this->geo == 1) {
        this->mode = 2;
      }
    } else if (this->fe2 == 1) {
      if (this->geo == 0) {
        this->mode = 101;
      } else if (this->geo == 1) {
        this->mode = 102;
      }
    }
  } else if (thermal == 1) {
    if (this->fe2 == 0) {
      if (this->geo == 0) {
        this->mode = 11;
      } else if (this->geo == 1) {
        this->mode = 12;
      }
    } else if (this->fe2 == 1) {
      if (this->geo == 0) {
        this->mode = 111;
      } else if (this->geo == 1) {
        this->mode = 112;
      }
    }
  }
  if (this->fe2 == 0) {
    this->E = list.getPrecVal("E");
    this->G = list.getPrecVal("G");
    this->A = list.getPrecVal("A");
    this->Ix = list.getPrecVal("Ix");
    this->Iy = list.getPrecVal("Iy");
    this->Iz = list.getPrecVal("Iz");
    this->ky = list.getPrecVal("ky");
    this->kz = list.getPrecVal("kz");
    this->alpha = list.getPrecVal("alpha");
    this->kappa = list.getPrecVal("kappa");
    this->RVECmat = list.getIndexVal("RVECmat");
  }

  this->meshIdDisp = list.getIndexVal("meshidDisp");
  this->meshIdRot = list.getIndexVal("meshidRot");
  this->meshIdTemp = list.getIndexVal("meshidTemp");

  this->disporder = list.getIndexVal("dispOrder");
  this->rotorder = list.getIndexVal("rotOrder");
  this->temporder = list.getIndexVal("tempOrder");

  if (this->fe2 == 0) {
    auto &Logger = pointers.getSPDLogger();

    Logger.info(
        "\n{:-<100}\n"
        "*   Specified parameters for element formulation 103 Timoshenko 3D:\n"
        "       Youngs's modulus E:         {:>12.4e}"
        "       Shear modulus G:            {:>12.4e}"
        "       Cross-section area A:       {:>12.4e}"
        "       Torsional stiffness It:     {:>12.4e}"
        "       Second moment of area Iy:   {:>12.4e}"
        "       Second moment of area Iz:   {:>12.4e}"
        "       Shear correction factor ky: {:>12.4e}"
        "       Shear correction factor kz: {:>12.4e}",
        "", this->E, this->G, this->A, this->Ix, this->Iy, this->Iz, this->ky,
        this->kz);

  } else {
    auto &Logger = pointers.getSPDLogger();

    Logger.info("Specified parameters for element formulation 103 Timoshenko "
                "3D:\n   Using the FE2 scheme.");
  }
  this->messageUnprocessed(pointers, list, "EL103_Timoshenko3D");
}

void EL103_Timoshenko3D::setDegreesOfFreedom(PointerCollection &pointers,
                                             FiniteElement::Edge &elem) {

  elem.setH1Shapes(pointers, meshIdDisp, disporder);
  elem.setH1Shapes(pointers, meshIdRot, rotorder);
  if (this->thermal == 1) {
    elem.setH1Shapes(pointers, meshIdTemp, temporder);
  }
}

void EL103_Timoshenko3D::AdditionalOperations(
    PointerCollection &pointers, FiniteElement::Edge &elem) {
  if (this->fe2 == 1) {
    auto GP = this->getIntegrationPoints(pointers, &elem);
    for (auto i : GP) {
      elem.getMaterialFormulation(pointers)->initRVE(pointers, i);
    }
  } else {
    if (this->RVECmat == 1 && !m_RVECmatInit) {
      m_RVECmatInit = true;
      auto GP = this->getIntegrationPoints(pointers, &elem);
      GP.setOrder(0);

      elem.getMaterialFormulation(pointers)->initRVE(pointers,
                                                      GP.getIntegrationPoint());

      Materials::MaterialTransferData materialData;
      if (thermal == 1) {
        materialData.strains.resize(12);
      } else {
        materialData.strains.resize(6);
      }

      materialData.strains.setZero();
      elem.getMaterialFormulation(pointers)->getMaterialData(
          pointers, materialData, GP.getIntegrationPoint());
      this->m_CMat = materialData.materialTangent;
      pointers.getSPDLogger().info(
          "Element 103 Timoshenko 3D, computed Material parameters:\n{}",
          m_CMat);
    } else if (this->RVECmat == 0) {
      if (thermal == 1) {
        this->m_CMat = this->getLinearCMatThermal();
      } else {
        this->m_CMat = this->getLinearCMat();
      }
    }
  }
}

void EL103_Timoshenko3D::updateRVEHistory(
    PointerCollection &pointers, FiniteElement::GenericFiniteElement *elem) {
  if (this->fe2) {
    auto GP = this->getIntegrationPoints(pointers, elem);
    for (auto i : GP) {
      elem->getMaterialFormulation(pointers)->updateRVEHistory(pointers, i);
    }
  }
}

auto EL103_Timoshenko3D::getDofs(PointerCollection &pointers,
                                 FiniteElement::GenericFiniteElement *elem)
    -> std::vector<DegreeOfFreedom *> {

  std::vector<DegreeOfFreedom *> Dofs, dispDofs, rotDofs;
  elem->getH1Dofs(pointers, dispDofs, this->meshIdDisp, this->disporder);
  elem->getH1Dofs(pointers, rotDofs, this->meshIdRot, this->rotorder);
  Dofs.clear();
  Dofs.insert(Dofs.end(), dispDofs.begin(), dispDofs.end());
  Dofs.insert(Dofs.end(), rotDofs.begin(), rotDofs.end());
  if (this->thermal == 1) {
    dispDofs.clear();
    elem->getH1Dofs(pointers, dispDofs, meshIdTemp, temporder);
    Dofs.insert(Dofs.end(), dispDofs.begin(), dispDofs.end());
  }
  return Dofs;
}

void EL103_Timoshenko3D::setTangentResidual(
    PointerCollection &pointers, FiniteElement::Edge &elem,
    Eigen::Matrix<prec, Eigen::Dynamic, Eigen::Dynamic> &stiffness,
    Eigen::Matrix<prec, Eigen::Dynamic, 1> &residual,
    std::vector<DegreeOfFreedom *> &Dofs) {

  switch (this->mode) {
  case 1: {
    this->setTangentResidualLinear(pointers, elem, stiffness, residual, Dofs);
  } break;
  case 101: {
    this->setTangentResidualLinearFE2(pointers, elem, stiffness, residual,
                                      Dofs);
  } break;
  case 2: {
    this->setTangentResidualNonlinear(pointers, elem, stiffness, residual,
                                      Dofs);
  } break;
  case 11: {
    this->setTangentResidualLinearThermal(pointers, elem, stiffness, residual,
                                          Dofs);
  } break;
  case 111: {
    this->setTangentResidualLinearThermalFE2(pointers, elem, stiffness,
                                             residual, Dofs);
  } break;
  default: {
    throw std::runtime_error("\nSelected mode for EL102 not implemented!");
  }
  }
}

void EL103_Timoshenko3D::setTangentResidualLinear(
    PointerCollection &pointers, FiniteElement::Edge &elem,
    Eigen::Matrix<prec, Eigen::Dynamic, Eigen::Dynamic> &stiffness,
    Eigen::Matrix<prec, Eigen::Dynamic, 1> &residual,
    std::vector<DegreeOfFreedom *> &Dofs) {

  Dofs = this->getDofs(pointers, &elem);
  auto GP = elem.getIntegrationPoints(pointers);

  indexType maxOrder = std::max({disporder, rotorder});

  GP.setOrder(maxOrder * 2);
  Types::MatrixXX<prec> B;
  B.resize(6, Dofs.size());
  B.setZero();
  stiffness.resize(Dofs.size(), Dofs.size());
  stiffness.setZero();
  residual.resize(Dofs.size());
  residual.setZero();
  auto solution = elem.getSolution(pointers, Dofs);
  for (auto &gp : GP) {
    auto jaco = elem.getJacobian(pointers, gp);
    auto shapesDisp = elem.getH1Shapes(pointers, disporder, jaco, gp);
    auto shapesRot = elem.getH1Shapes(pointers, rotorder, jaco, gp);
    indexType pos = 0;
    for (auto i = 0; i < shapesDisp.shapes.rows(); ++i) {
      B(0, pos) = shapesDisp.shapeDeriv(0, i);
      B(1, pos + 1) = shapesDisp.shapeDeriv(0, i);
      B(2, pos + 2) = shapesDisp.shapeDeriv(0, i);
      pos += 3;
    }
    for (auto i = 0; i < shapesRot.shapes.rows(); ++i) {
      B(1, pos + 2) = shapesRot.shapes(i);
      B(2, pos + 1) = shapesRot.shapes(i);
      B(3, pos + 0) = shapesRot.shapeDeriv(0, i);
      B(4, pos + 1) = shapesRot.shapeDeriv(0, i);
      B(5, pos + 2) = shapesRot.shapeDeriv(0, i);
      pos += 3;
    }
    stiffness += B.transpose() * m_CMat * B * gp.weight * jaco;
  }
  residual = stiffness * solution;
}

void EL103_Timoshenko3D::setTangentResidualLinearFE2(
    PointerCollection &pointers, FiniteElement::Edge &elem,
    Eigen::Matrix<prec, Eigen::Dynamic, Eigen::Dynamic> &stiffness,
    Eigen::Matrix<prec, Eigen::Dynamic, 1> &residual,
    std::vector<DegreeOfFreedom *> &Dofs) {

  Dofs = this->getDofs(pointers, &elem);
  auto GP = elem.getIntegrationPoints(pointers);

  indexType maxOrder = std::max({disporder, rotorder});

  GP.setOrder(maxOrder * 2);
  Types::MatrixXX<prec> B;
  B.resize(6, Dofs.size());
  B.setZero();
  stiffness.resize(Dofs.size(), Dofs.size());
  stiffness.setZero();
  residual.resize(Dofs.size());
  residual.setZero();
  auto solution = elem.getSolution(pointers, Dofs);

  Materials::MaterialTransferData materialData;
  auto Hist = elem.getHistoryDataIterator(pointers);
  for (auto &gp : GP) {
    auto jaco = elem.getJacobian(pointers, gp);
    auto shapesDisp = elem.getH1Shapes(pointers, disporder, jaco, gp);
    auto shapesRot = elem.getH1Shapes(pointers, rotorder, jaco, gp);
    indexType pos = 0;
    for (auto i = 0; i < shapesDisp.shapes.rows(); ++i) {
      B(0, pos) = shapesDisp.shapeDeriv(0, i);
      B(1, pos + 1) = shapesDisp.shapeDeriv(0, i);
      B(2, pos + 2) = shapesDisp.shapeDeriv(0, i);
      pos += 3;
    }
    for (auto i = 0; i < shapesRot.shapes.rows(); ++i) {
      B(1, pos + 2) = shapesRot.shapes(i);
      B(2, pos + 1) = shapesRot.shapes(i);
      B(3, pos + 0) = shapesRot.shapeDeriv(0, i);
      B(4, pos + 1) = shapesRot.shapeDeriv(0, i);
      B(5, pos + 2) = shapesRot.shapeDeriv(0, i);
      pos += 3;
    }
    materialData.strains = B * solution;
    elem.getMaterialFormulation(pointers)->getMaterialData(pointers,
                                                           materialData, gp);
    // std::cout << "\nStrains: " << materialData.strains.transpose();
    // std::cout << "\nStresses: " << materialData.stresses.transpose();
    // std::cout << "\nMaterial tangent: \n"
    //           << materialData.materialTangent << "\n"
    //           << std::endl;
    stiffness +=
        B.transpose() * materialData.materialTangent * B * gp.weight * jaco;
    residual += B.transpose() * materialData.stresses * gp.weight * jaco;

    Hist.next();
  }
}

void EL103_Timoshenko3D::setTangentResidualLinearThermal(
    PointerCollection &pointers, FiniteElement::Edge &elem,
    Eigen::Matrix<prec, Eigen::Dynamic, Eigen::Dynamic> &stiffness,
    Eigen::Matrix<prec, Eigen::Dynamic, 1> &residual,
    std::vector<DegreeOfFreedom *> &Dofs) {

  Dofs = this->getDofs(pointers, &elem);
  auto GP = elem.getIntegrationPoints(pointers);

  indexType maxOrder = std::max({disporder, rotorder, temporder});

  GP.setOrder(maxOrder * 2);
  Types::MatrixXX<prec> B;
  B.resize(12, Dofs.size());
  B.setZero();
  stiffness.resize(Dofs.size(), Dofs.size());
  stiffness.setZero();
  residual.resize(Dofs.size());
  residual.setZero();
  auto solution = elem.getSolution(pointers, Dofs);
  for (auto &gp : GP) {
    auto jaco = elem.getJacobian(pointers, gp);
    auto shapesDisp = elem.getH1Shapes(pointers, disporder, jaco, gp);
    auto shapesRot = elem.getH1Shapes(pointers, rotorder, jaco, gp);
    auto shapesThermal = elem.getH1Shapes(pointers, temporder, jaco, gp);
    indexType pos = 0;
    for (auto i = 0; i < shapesDisp.shapes.rows(); ++i) {
      B(0, pos) = shapesDisp.shapeDeriv(0, i);
      B(1, pos + 1) = shapesDisp.shapeDeriv(0, i);
      B(2, pos + 2) = shapesDisp.shapeDeriv(0, i);
      pos += 3;
    }
    for (auto i = 0; i < shapesRot.shapes.rows(); ++i) {
      B(1, pos + 2) = shapesRot.shapes(i);
      B(2, pos + 1) = shapesRot.shapes(i);
      B(3, pos + 0) = shapesRot.shapeDeriv(0, i);
      B(4, pos + 1) = shapesRot.shapeDeriv(0, i);
      B(5, pos + 2) = shapesRot.shapeDeriv(0, i);
      pos += 3;
    }
    for (auto i = 0; i < shapesThermal.shapes.rows(); ++i) {
      B(6, pos) = shapesDisp.shapes(i);
      B(7, pos + 1) = shapesDisp.shapes(i);
      B(8, pos + 2) = shapesDisp.shapes(i);
      B(9, pos) = shapesDisp.shapeDeriv(0, i);
      B(10, pos + 1) = shapesDisp.shapeDeriv(0, i);
      B(11, pos + 2) = shapesDisp.shapeDeriv(0, i);
      pos += 3;
    }
    stiffness += B.transpose() * m_CMat * B * gp.weight * jaco;
  }
  residual = stiffness * solution;
}

void EL103_Timoshenko3D::setTangentResidualLinearThermalFE2(
    PointerCollection &pointers, FiniteElement::Edge &elem,
    Eigen::Matrix<prec, Eigen::Dynamic, Eigen::Dynamic> &stiffness,
    Eigen::Matrix<prec, Eigen::Dynamic, 1> &residual,
    std::vector<DegreeOfFreedom *> &Dofs) {

  Dofs = this->getDofs(pointers, &elem);
  auto GP = elem.getIntegrationPoints(pointers);

  indexType maxOrder = std::max({disporder, rotorder, temporder});

  GP.setOrder(maxOrder * 2);
  Types::MatrixXX<prec> B;
  B.resize(12, Dofs.size());
  B.setZero();
  stiffness.resize(Dofs.size(), Dofs.size());
  stiffness.setZero();
  residual.resize(Dofs.size());
  residual.setZero();
  auto solution = elem.getSolution(pointers, Dofs);

  auto Hist = elem.getHistoryDataIterator(pointers);
  Materials::MaterialTransferData materialData;
  materialData.strains.resize(12);
  materialData.historyData = &Hist;
  for (auto &gp : GP) {
    auto jaco = elem.getJacobian(pointers, gp);
    auto shapesDisp = elem.getH1Shapes(pointers, disporder, jaco, gp);
    auto shapesRot = elem.getH1Shapes(pointers, rotorder, jaco, gp);
    auto shapesThermal = elem.getH1Shapes(pointers, temporder, jaco, gp);
    indexType pos = 0;
    for (auto i = 0; i < shapesDisp.shapes.rows(); ++i) {
      B(0, pos) = shapesDisp.shapeDeriv(0, i);
      B(1, pos + 1) = shapesDisp.shapeDeriv(0, i);
      B(2, pos + 2) = shapesDisp.shapeDeriv(0, i);
      pos += 3;
    }
    for (auto i = 0; i < shapesRot.shapes.rows(); ++i) {
      B(1, pos + 2) = shapesRot.shapes(i);
      B(2, pos + 1) = shapesRot.shapes(i);
      B(3, pos + 0) = shapesRot.shapeDeriv(0, i);
      B(4, pos + 1) = shapesRot.shapeDeriv(0, i);
      B(5, pos + 2) = shapesRot.shapeDeriv(0, i);
      pos += 3;
    }
    for (auto i = 0; i < shapesThermal.shapes.rows(); ++i) {
      B(6, pos) = shapesDisp.shapes(i);
      B(7, pos + 1) = shapesDisp.shapes(i);
      B(8, pos + 2) = shapesDisp.shapes(i);
      B(9, pos) = shapesDisp.shapeDeriv(0, i);
      B(10, pos + 1) = shapesDisp.shapeDeriv(0, i);
      B(11, pos + 2) = shapesDisp.shapeDeriv(0, i);
      pos += 3;
    }

    materialData.strains = B * solution;
    elem.getMaterialFormulation(pointers)->getMaterialData(pointers,
                                                           materialData, gp);

    stiffness +=
        B.transpose() * materialData.materialTangent * B * gp.weight * jaco;
    residual += B.transpose() * materialData.stresses * gp.weight * jaco;
    Hist.next();
  }
}

void EL103_Timoshenko3D::setTangentResidualNonlinear(
    PointerCollection &pointers, FiniteElement::Edge &elem,
    Eigen::Matrix<prec, Eigen::Dynamic, Eigen::Dynamic> &stiffness,
    Eigen::Matrix<prec, Eigen::Dynamic, 1> &residual,
    std::vector<DegreeOfFreedom *> &Dofs) {}

auto EL103_Timoshenko3D::getLinearCMat() -> Types::Matrix66<prec> {
  Types::Matrix66<prec> C;
  C.setZero();

  C(0, 0) = E * A;
  C(1, 1) = G * A * ky;
  C(2, 2) = G * A * kz;
  C(3, 3) = G * Ix;
  C(4, 4) = E * Iy;
  C(5, 5) = E * Iz;

  return C;
}

auto EL103_Timoshenko3D::getLinearCMatThermal() -> Types::MatrixXX<prec> {
  Types::MatrixXX<prec> C;
  C.resize(12, 12);
  C.setZero();

  C(0, 0) = E * A;
  C(1, 1) = G * A * ky;
  C(2, 2) = G * A * kz;
  C(3, 3) = G * Ix;
  C(4, 4) = E * Iy;
  C(5, 5) = E * Iz;
  C(0, 6) = -E * A * alpha;
  C(4, 7) = -E * Iy * alpha;
  C(5, 8) = -E * Iz * alpha;

  C(7, 7) = -kappa * A;
  C(8, 8) = -kappa * A;
  C(9, 9) = -kappa * A;
  C(10, 10) = -kappa * Iy;
  C(11, 11) = -kappa * Iz;

  return C;
}

auto EL103_Timoshenko3D::getIntegrationPoints(
    PointerCollection &pointers, FiniteElement::GenericFiniteElement *elem)
    -> IntegrationPoints {

  auto GP = elem->getIntegrationPoints(pointers);
  indexType mOrd = std::max(rotorder, disporder);
  mOrd = std::max(mOrd, temporder);
  GP.setOrder(mOrd * 2);
  return GP;
}

void EL103_Timoshenko3D::setMass(
    PointerCollection &pointers, FiniteElement::GenericFiniteElement *elem,
    Eigen::Matrix<prec, Eigen::Dynamic, Eigen::Dynamic> &stiffness,
    Eigen::Matrix<prec, Eigen::Dynamic, 1> &residual,
    std::vector<DegreeOfFreedom *> &Dofs) {}

void EL103_Timoshenko3D::getElementsLocalNodalReactions(
    PointerCollection &ptrCol, FiniteElement::GenericFiniteElement *elem,
    std::map<indexType, std::vector<prec>> &vReacs) {}

void EL103_Timoshenko3D::toParaviewAdaper(PointerCollection &pointers,
                                          FiniteElement::Edge &elem,
                                          vtkPlotInterface &paraviewAdapter,
                                          ParaviewSwitch control) {

  int matNum = static_cast<int>(elem.getMaterial()->getNumber());
  switch (control) {
  case HierAMuS::ParaviewSwitch::Mesh: {
    elem.geometryToParaview(pointers, paraviewAdapter, 0, matNum);
    break;
  }
  case HierAMuS::ParaviewSwitch::Solution: {
    indexType maxOrder = std::max({disporder, rotorder, temporder});
    elem.H1SolutionToParaview(pointers, paraviewAdapter, 0, matNum,
                               this->meshIdDisp, maxOrder,
                               paraviewNames::DisplacementName());

    elem.H1SolutionToParaview(pointers, paraviewAdapter, 0, matNum,
                               this->meshIdRot, maxOrder,
                               paraviewNames::RotationName());
    if (this->thermal == 1) {
      elem.H1SolutionToParaview(pointers, paraviewAdapter, 0, matNum,
                                 this->meshIdTemp, maxOrder,
                                 paraviewNames::TemperatureName());
    }
    break;
  }
  case ParaviewSwitch::Weights: {
    elem.computeWeightsParaview(pointers, paraviewAdapter, 0, matNum);
  } break;
  default: {
    break;
  }
  }
}

auto EL103_Timoshenko3D::getHistoryDataStructure()
    -> const HistoryDataStructure & {
  return m_HistoryDataStructure;
}

auto EL103_Timoshenko3D::getNumberOfIntergrationPoints(
    PointerCollection &pointers, FiniteElement::GenericFiniteElement *elem)
    -> indexType {
  auto GP = elem->getIntegrationPoints(pointers);
  indexType order = std::max(this->disporder, this->rotorder);
  GP.setOrder(order * 2);
  return GP.getTotalGP();
}

Types::Matrix33<prec> EL103_Timoshenko3D::getMaterialMatrix() const {
  Types::Matrix33<prec> mat;
  return mat;
}

const HistoryDataStructure
    EL103_Timoshenko3D::m_HistoryDataStructure({{4, 2}, {8, 1}}, {});

} // namespace HierAMuS::Elementformulations
