// Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause



#include "plot/vtkplotClassBase.h"
#include "shapefunctions/IntegrationsPoints/IntegrationPoints.h"


#include <elementFormulations/EL206_Plate.h>

#include <control/HandlingStructs.h>
#include <control/OutputHandler.h>
#include "control/ParameterList.h"
#include <pointercollection/pointercollection.h>

#include <finiteElements/Face.h>

#include <materials/GenericMaterialFormulation.h>

#include <geometry/VertexData.h>

#include <elementFormulations/GenericElementFormulation.h>
#include <solver/GenericSolutionState.h>

#include <types/MatrixTypes.h>

#include <Eigen/Eigenvalues>
#include <vector>

#include <vtkCellType.h>

namespace HierAMuS::Elementformulations {

EL206_Plate::EL206_Plate(PointerCollection *ptrCol)
    : GenericElementFormulationInterface(ptrCol) {}

EL206_Plate::~EL206_Plate() = default;

void EL206_Plate::readData(PointerCollection &pointers, ParameterList &list) {

  this->meshIdDisp = list.getIndexVal("meshiddisp");
  this->meshIdRot = list.getIndexVal("meshidrot");
  this->dispOrder = list.getIndexVal("orderdisp");
  this->rotOrder = list.getIndexVal("orderrot");
  this->mode = list.getIndexVal("mode");
  this->EI = list.getPrecVal("EI");
  this->GA = list.getPrecVal("GA");
  this->GI = list.getPrecVal("GI");

  auto &Log = pointers.getSPDLogger();

  Log.info("\n{:-<100}\n"
                "*   Element 206, specified Options\n"
                "    Mesh id for displacement nodes:         {:>12}\n"
                "    Shape function order for displacements: {:>12}\n"
                "{:-<100}\n",
                "",
                this->meshIdDisp,
                this->dispOrder,
                "");

  this->messageUnprocessed(pointers, list, "EL206_Plate");
}

void EL206_Plate::setDegreesOfFreedom(
  PointerCollection& pointers, FiniteElement::Face &elem) {

  switch (this->mode) {
  case 1: {
    elem.setH1Shapes(pointers, this->meshIdDisp, this->dispOrder);
    elem.setH1Shapes(pointers, this->meshIdRot, this->rotOrder);
  } break;
  case 2: {
    elem.setSpecialPlateShapes(pointers, this->meshIdRot, this->rotOrder);
  } break;
  }
}

void EL206_Plate::AdditionalOperations(
  PointerCollection& pointers, FiniteElement::Face &elem) {

  switch (this->mode) {
  case 1: {
    elem.setAllNodeBoundaryConditionMeshId(pointers, this->meshIdDisp, 0);
    elem.setAllNodeBoundaryConditionMeshId(pointers, this->meshIdDisp, 1);
    elem.setAllNodeBoundaryConditionMeshId(pointers, this->meshIdRot, 2);
  }
  case 2: {
    elem.setAllNodeBoundaryConditionMeshId(pointers, this->meshIdDisp, 0);
    elem.setAllNodeBoundaryConditionMeshId(pointers, this->meshIdDisp, 1);
    elem.setAllNodeBoundaryConditionMeshId(pointers, this->meshIdRot, 1);
    elem.setAllNodeBoundaryConditionMeshId(pointers, this->meshIdRot, 2);
  }
  }
}

auto EL206_Plate::getDofs(PointerCollection& pointers, FiniteElement::GenericFiniteElement *elem)
-> std::vector<DegreeOfFreedom *> {

  std::vector<DegreeOfFreedom *> Dofs;
  Dofs.clear();
  std::vector<DegreeOfFreedom *> tempDofs;
  elem->getH1Dofs(pointers, tempDofs, this->meshIdDisp, this->dispOrder);
  Dofs.insert(Dofs.end(), tempDofs.begin(), tempDofs.end());
  tempDofs.clear();
  elem->getH1Dofs(pointers, tempDofs, this->meshIdRot, this->rotOrder);
  Dofs.insert(Dofs.end(), tempDofs.begin(), tempDofs.end());
  return Dofs;
}

void EL206_Plate::setTangentResidual(PointerCollection& pointers,
                                     FiniteElement::Face &elem,
                                     Types::MatrixXX<prec> &stiffness,
                                     Types::VectorX<prec> &residual, std::vector<DegreeOfFreedom *> &Dofs) {

  switch (this->mode) {
  case 1: // Displacement based formulation
    this->setTangentResidualDispFormulation(pointers, &elem, stiffness, residual,
                                            Dofs);
    break;
  case 2: // Hu-Washizu formulation
    this->setTangentResidualHuWashizuFormulation(pointers, &elem, stiffness,
                                                   residual, Dofs);
    break;
  default:
    throw std::runtime_error("Element formulation 201 called with wrong mode!");
  }
}




auto EL206_Plate::getNumberOfIntergrationPoints(
    PointerCollection &pointers, FiniteElement::GenericFiniteElement *elem)
    -> indexType {
  auto GP = elem->getIntegrationPoints(pointers);
  GP.setOrder(2);
  return GP.getTotalGP();
}

void EL206_Plate::toParaviewAdaper(PointerCollection &pointers,
                                   FiniteElement::Face &elem,
                                   vtkPlotInterface &paraviewAdapter,
                                   ParaviewSwitch control) {

  // int matNum = static_cast<int>(elem->getMaterial()->getNumber());
  // switch (control) {
  //   case ParaviewSwitch::Mesh: {
  //   elem->geometryToParaview(pointers, paraviewAdapter, 0, matNum);

  //   }break;
  //   case ParaviewSwitch::Solution: {
  //     elem->H1SolutionToParaview(*this->ptrCol, paraviewAdapter, 0, matNum,
  //     this->meshIdDisp, this->intOrderDisp,
  //     paraviewNames::DisplacementName());
  //   }break;
  //   case ParaviewSwitch::Weights: {
  //     elem->computeWeightsParaview(pointers, paraviewAdapter, 0, matNum);
  //   } break;
  //  }
  // switch (elem->getType()) {
  // case FiniteElement::Elementtypes::Quadrilateral: {
  //   this->toParaviewQuadrilateral(pointers, elem, paraviewAdapter, control);
  // } break;
  // default:
  //   throw std::runtime_error("Paraview not implemented in this "
  //                            "elementformulation for the given element
  //                            type.");
  // }
}

void EL206_Plate::toParaviewQuadrilateral(
    PointerCollection &pointers, FiniteElement::GenericFiniteElement *elem,
    vtkPlotInterface &paraviewAdapter, ParaviewSwitch control) {}

void EL206_Plate::setTangentResidualDispFormulation(
  PointerCollection& pointers,
  FiniteElement::Face *elem,
  Eigen::Matrix<prec, Eigen::Dynamic, Eigen::Dynamic> &stiffness,
  Eigen::Matrix<prec, Eigen::Dynamic, 1> &residual, std::vector<DegreeOfFreedom *> &Dofs) {

  Dofs.clear();
  std::vector<DegreeOfFreedom *> tempDofs;
  elem->getH1Dofs(pointers, tempDofs, this->meshIdDisp, this->dispOrder);
  Dofs.insert(Dofs.end(), tempDofs.begin(), tempDofs.end());
  tempDofs.clear();
  elem->getH1Dofs(pointers, tempDofs, this->meshIdRot, this->rotOrder);
  Dofs.insert(Dofs.end(), tempDofs.begin(), tempDofs.end());

  indexType nDofs = Dofs.size();

  stiffness.resize(nDofs, nDofs);
  residual.resize(nDofs);
  stiffness.setZero();
  residual.setZero();

  auto GP = elem->getIntegrationPoints(pointers);
  GP.setOrder(2);
  for (auto i : GP) {
    auto jacob = elem->getJacobian(pointers, i);
    auto shapesDisp = elem->getH1Shapes(pointers, this->dispOrder, jacob, i);
    auto shapesRot = elem->getH1Shapes(pointers, this->rotOrder, jacob, i);

    auto Bmatrix = this->getBMatrix(shapesDisp, shapesRot);
    std::cout << i.xi << std::endl;
  }
}


Types::Matrix33<prec> EL206_Plate::getMaterialMatrix() {
  Types::Matrix3X<prec> material;
  material.resize(5,5);
  material.setZero();
  return {};
}

Types::Matrix3X<prec> EL206_Plate::getBMatrix(Geometry::H1Shapes &shapesDisp,
                                              Geometry::H1Shapes &shapesRot) {

  Types::Matrix3X<prec> Bmat;
  indexType ndispShapes = shapesDisp.shapes.size();
  indexType nrotShapes = shapesRot.shapes.size();

  Bmat.resize(5, (ndispShapes + nrotShapes) * 3);
  Bmat.setZero();

  return Bmat;
}

void EL206_Plate::setTangentResidualHuWashizuFormulation(
  PointerCollection& pointers,
  FiniteElement::Face *elem,
  Eigen::Matrix<prec, Eigen::Dynamic, Eigen::Dynamic> &stiffness,
  Eigen::Matrix<prec, Eigen::Dynamic, 1> &residual, std::vector<DegreeOfFreedom *> &Dofs) {}

Types::Matrix3X<prec> EL206_Plate::getLocalStrainStressInterpolation(
    FiniteElement::GenericFiniteElement *elem, prec xi, prec eta) {

  return {};
}

} // namespace HierAMuS
