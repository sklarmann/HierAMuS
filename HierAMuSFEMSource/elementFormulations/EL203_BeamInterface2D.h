// Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause



#pragma once

#include <forwarddeclaration.h>

#include <types/MatrixTypes.h>

#include <materials/GenericMaterialFormulation.h>

#include <Eigen/Dense>
#include <elementFormulations/GenericElementFormulation.h>
#include <vector>

namespace HierAMuS::Elementformulations {

class EL203_BeamInterface2D : public GenericElementFormulation {
public:
  explicit EL203_BeamInterface2D(PointerCollection *ptrCol);
  ~EL203_BeamInterface2D() override;
  void readData(PointerCollection &pointers, ParameterList &list) override;
  void setDegreesOfFreedom(PointerCollection& pointers, FiniteElement::GenericFiniteElement *elem) override;
  void AdditionalOperations(PointerCollection& pointers, FiniteElement::GenericFiniteElement *elem) override;
  auto getDofs(PointerCollection& pointers, FiniteElement::GenericFiniteElement *elem)
  -> std::vector<DegreeOfFreedom *> override;
  void setTangentResidual(
    PointerCollection& pointers,
    FiniteElement::GenericFiniteElement *elem,
    Eigen::Matrix<prec, Eigen::Dynamic, Eigen::Dynamic> &stiffness,
    Eigen::Matrix<prec, Eigen::Dynamic, 1> &residual, std::vector<DegreeOfFreedom *> &Dofs) override;

  auto getNumberOfIntergrationPoints(PointerCollection &pointers,
                                     FiniteElement::GenericFiniteElement *elem)
      -> indexType override;

  void toParaviewAdaper(PointerCollection &pointers,
                        FiniteElement::GenericFiniteElement *elem,
                        vtkPlotInterface &paraviewAdapter,
                        ParaviewSwitch control) override;

private:
  // Mode 1
  void setDegreesOfFreedomV1(PointerCollection& pointers, FiniteElement::beamInterfaceElement2D *elem);
  void AdditionalOperationsV1(PointerCollection& pointers, FiniteElement::beamInterfaceElement2D *elem);
  void setTangentResidualV1(
    PointerCollection& pointers,
    FiniteElement::beamInterfaceElement2D *elem,
    Eigen::Matrix<prec, Eigen::Dynamic, Eigen::Dynamic> &stiffness,
    Eigen::Matrix<prec, Eigen::Dynamic, 1> &residual, std::vector<DegreeOfFreedom *> &Dofs);
  // Mode 2
  void setDegreesOfFreedomV2(PointerCollection& pointers, FiniteElement::beamInterfaceElement2D *elem);
  void AdditionalOperationsV2(PointerCollection& pointers, FiniteElement::beamInterfaceElement2D *elem);
  void setTangentResidualV2(
    PointerCollection& pointers,
    FiniteElement::beamInterfaceElement2D *elem,
    Eigen::Matrix<prec, Eigen::Dynamic, Eigen::Dynamic> &stiffness,
    Eigen::Matrix<prec, Eigen::Dynamic, 1> &residual, std::vector<DegreeOfFreedom *> &Dofs);

  // Mode 3 - Mixed Formulation
  void setDegreesOfFreedomV3(PointerCollection& pointers, FiniteElement::beamInterfaceElement2D *elem);
  void AdditionalOperationsV3(PointerCollection& pointers, FiniteElement::beamInterfaceElement2D *elem);
  void setTangentResidualV3(
    PointerCollection& pointers,
    FiniteElement::beamInterfaceElement2D *elem,
    Eigen::Matrix<prec, Eigen::Dynamic, Eigen::Dynamic> &stiffness,
    Eigen::Matrix<prec, Eigen::Dynamic, 1> &residual, std::vector<DegreeOfFreedom *> &Dofs);

  // Mode 4 - Mixed Formulation, local stress approximation
  void setDegreesOfFreedomV4(PointerCollection& pointers, FiniteElement::beamInterfaceElement2D *elem);
  void AdditionalOperationsV4(PointerCollection& pointers, FiniteElement::beamInterfaceElement2D *elem);
  void setTangentResidualV4(
    PointerCollection& pointers,
    FiniteElement::beamInterfaceElement2D *elem,
    Eigen::Matrix<prec, Eigen::Dynamic, Eigen::Dynamic> &stiffness,
    Eigen::Matrix<prec, Eigen::Dynamic, 1> &residual, std::vector<DegreeOfFreedom *> &Dofs);

  // Mode 5 - Mixed Formulation, local stress approximation
  void setDegreesOfFreedomV5(PointerCollection& pointers, FiniteElement::beamInterfaceElement2D &elem);
  void AdditionalOperationsV5(PointerCollection& pointers, FiniteElement::beamInterfaceElement2D &elem);
  void setTangentResidualV5(
    PointerCollection& pointers,
    FiniteElement::beamInterfaceElement2D *elem,
    Eigen::Matrix<prec, Eigen::Dynamic, Eigen::Dynamic> &stiffness,
    Eigen::Matrix<prec, Eigen::Dynamic, 1> &residual, std::vector<DegreeOfFreedom *> &Dofs);

  // Mode 6 - Mixed Formulation, local stress approximation
  void setDegreesOfFreedomV6(PointerCollection& pointers, FiniteElement::beamInterfaceElement2D &elem);
  void AdditionalOperationsV6(PointerCollection& pointers, FiniteElement::beamInterfaceElement2D &elem);
  void setTangentResidualV6(
    PointerCollection& pointers,
    FiniteElement::beamInterfaceElement2D *elem,
    Eigen::Matrix<prec, Eigen::Dynamic, Eigen::Dynamic> &stiffness,
    Eigen::Matrix<prec, Eigen::Dynamic, 1> &residual, std::vector<DegreeOfFreedom *> &Dofs);

  
  // Mode 7 - Mixed Formulation, local stress approximation
  void setDegreesOfFreedomV7(PointerCollection& pointers, FiniteElement::beamInterfaceElement2D &elem);
  void AdditionalOperationsV7(PointerCollection& pointers, FiniteElement::beamInterfaceElement2D &elem);
  void setTangentResidualV7(
    PointerCollection& pointers,
    FiniteElement::beamInterfaceElement2D *elem,
    Eigen::Matrix<prec, Eigen::Dynamic, Eigen::Dynamic> &stiffness,
    Eigen::Matrix<prec, Eigen::Dynamic, 1> &residual, std::vector<DegreeOfFreedom *> &Dofs);

  void computeStressesAndStrains(FiniteElement::beamInterfaceElement2D &elem,
                                 indexType localEdgeNumber, prec xi, prec eta);

  indexType meshIdDisp, intOrder;
  indexType meshIdWarp, meshIdRot;

  indexType mode;

  std::vector<indexType> edgeList;
  std::vector<indexType> materialList;
  indexType beamVertexNumber;

  std::map<indexType, indexType> edgeMaterialMap;
  std::map<indexType, indexType> plotVertexMap;
  std::map<indexType, indexType> plotVertexMapOpposite;

  std::vector<Materials::GenericMaterialFormulation *> Materials;
};

} // namespace HierAMuS
