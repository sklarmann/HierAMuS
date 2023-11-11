// Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include "elementFormulations/GenericElementFormulationInterface.h"
#include "shapefunctions/IntegrationsPoints/IntegrationPoints.h"
#include <types/MatrixTypes.h>

#include <Eigen/Dense>
#include <vector>

namespace HierAMuS::FiniteElement {
class beamInterfaceElement2D;
}

namespace HierAMuS::Materials {
class GenericMaterialFormulation;
}

namespace HierAMuS::Elementformulations {

class EL203_BeamInterface2D : public GenericElementFormulationInterface<FiniteElement::beamInterfaceElement2D> {
public:
  explicit EL203_BeamInterface2D(PointerCollection *ptrCol);
  ~EL203_BeamInterface2D() override;
  void readData(PointerCollection &pointers, ParameterList &list) override;
  void setDegreesOfFreedom(PointerCollection &pointers,
                           FiniteElement::beamInterfaceElement2D &elem) override;
  void AdditionalOperations(PointerCollection &pointers,
                            FiniteElement::beamInterfaceElement2D &elem) override;
  auto getDofs(PointerCollection &pointers,
               FiniteElement::GenericFiniteElement *elem)
      -> std::vector<DegreeOfFreedom *> override;
  void setTangentResidual(
      PointerCollection &pointers, FiniteElement::beamInterfaceElement2D &elem,
      Eigen::Matrix<prec, Eigen::Dynamic, Eigen::Dynamic> &stiffness,
      Eigen::Matrix<prec, Eigen::Dynamic, 1> &residual,
      std::vector<DegreeOfFreedom *> &Dofs) override;

  auto getNumberOfIntergrationPoints(PointerCollection &pointers,
                                     FiniteElement::GenericFiniteElement *elem)
      -> indexType override;
  auto getIntegrationPoints(PointerCollection &pointers,
                                     FiniteElement::GenericFiniteElement *elem)
      -> IntegrationPoints;

  void toParaviewAdaper(PointerCollection &pointers,
                        FiniteElement::beamInterfaceElement2D &elem,
                        vtkPlotInterface &paraviewAdapter,
                        ParaviewSwitch control) override;

  

  
private:
  // Mode 1
  void setDegreesOfFreedomV1(PointerCollection &pointers,
                             FiniteElement::beamInterfaceElement2D *elem);
  void AdditionalOperationsV1(PointerCollection &pointers,
                              FiniteElement::beamInterfaceElement2D *elem);
  void setTangentResidualV1(
      PointerCollection &pointers, FiniteElement::beamInterfaceElement2D *elem,
      Eigen::Matrix<prec, Eigen::Dynamic, Eigen::Dynamic> &stiffness,
      Eigen::Matrix<prec, Eigen::Dynamic, 1> &residual,
      std::vector<DegreeOfFreedom *> &Dofs);
  // Mode 2
  void setDegreesOfFreedomV2(PointerCollection &pointers,
                             FiniteElement::beamInterfaceElement2D *elem);
  void AdditionalOperationsV2(PointerCollection &pointers,
                              FiniteElement::beamInterfaceElement2D *elem);
  void setTangentResidualV2(
      PointerCollection &pointers, FiniteElement::beamInterfaceElement2D *elem,
      Eigen::Matrix<prec, Eigen::Dynamic, Eigen::Dynamic> &stiffness,
      Eigen::Matrix<prec, Eigen::Dynamic, 1> &residual,
      std::vector<DegreeOfFreedom *> &Dofs);

  // Mode 3 - Mixed Formulation
  void setDegreesOfFreedomV3(PointerCollection &pointers,
                             FiniteElement::beamInterfaceElement2D *elem);
  void AdditionalOperationsV3(PointerCollection &pointers,
                              FiniteElement::beamInterfaceElement2D *elem);
  void setTangentResidualV3(
      PointerCollection &pointers, FiniteElement::beamInterfaceElement2D *elem,
      Eigen::Matrix<prec, Eigen::Dynamic, Eigen::Dynamic> &stiffness,
      Eigen::Matrix<prec, Eigen::Dynamic, 1> &residual,
      std::vector<DegreeOfFreedom *> &Dofs);

  // Mode 4 - Mixed Formulation, local stress approximation
  void setDegreesOfFreedomV4(PointerCollection &pointers,
                             FiniteElement::beamInterfaceElement2D *elem);
  void AdditionalOperationsV4(PointerCollection &pointers,
                              FiniteElement::beamInterfaceElement2D *elem);
  void setTangentResidualV4(
      PointerCollection &pointers, FiniteElement::beamInterfaceElement2D *elem,
      Eigen::Matrix<prec, Eigen::Dynamic, Eigen::Dynamic> &stiffness,
      Eigen::Matrix<prec, Eigen::Dynamic, 1> &residual,
      std::vector<DegreeOfFreedom *> &Dofs);

  // Mode 5 - Mixed Formulation, local stress approximation
  void setDegreesOfFreedomV5(PointerCollection &pointers,
                             FiniteElement::beamInterfaceElement2D &elem);
  void AdditionalOperationsV5(PointerCollection &pointers,
                              FiniteElement::beamInterfaceElement2D &elem);
  void setTangentResidualV5(
      PointerCollection &pointers, FiniteElement::beamInterfaceElement2D *elem,
      Eigen::Matrix<prec, Eigen::Dynamic, Eigen::Dynamic> &stiffness,
      Eigen::Matrix<prec, Eigen::Dynamic, 1> &residual,
      std::vector<DegreeOfFreedom *> &Dofs);

  // Mode 6 - Mixed Formulation, local stress approximation
  void setDegreesOfFreedomV6(PointerCollection &pointers,
                             FiniteElement::beamInterfaceElement2D &elem);
  void AdditionalOperationsV6(PointerCollection &pointers,
                              FiniteElement::beamInterfaceElement2D &elem);
  void setTangentResidualV6(
      PointerCollection &pointers, FiniteElement::beamInterfaceElement2D *elem,
      Eigen::Matrix<prec, Eigen::Dynamic, Eigen::Dynamic> &stiffness,
      Eigen::Matrix<prec, Eigen::Dynamic, 1> &residual,
      std::vector<DegreeOfFreedom *> &Dofs);


  // Mode 7 - Mixed Formulation, local stress approximation
  void setDegreesOfFreedomV7(PointerCollection &pointers,
                             FiniteElement::beamInterfaceElement2D &elem);
  void AdditionalOperationsV7(PointerCollection &pointers,
                              FiniteElement::beamInterfaceElement2D &elem);
  void setTangentResidualV7(
      PointerCollection &pointers, FiniteElement::beamInterfaceElement2D *elem,
      Eigen::Matrix<prec, Eigen::Dynamic, Eigen::Dynamic> &stiffness,
      Eigen::Matrix<prec, Eigen::Dynamic, 1> &residual,
      std::vector<DegreeOfFreedom *> &Dofs);

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

  
  static constexpr indexType m_resEps_id = 0;
  static constexpr indexType m_resEpsC_id = 1;
  static constexpr indexType m_resSig_id = 2;

  static constexpr indexType m_Eps_id = 3;
  static constexpr indexType m_EpsC_id = 4;
  static constexpr indexType m_Sig_id = 5;

  static constexpr indexType m_G_id = 6;
  static constexpr indexType m_H_id = 7;
  static constexpr indexType m_J_id = 8;
  static constexpr indexType m_L_id = 9;
  static constexpr indexType m_M_id = 10;
};

} // namespace HierAMuS::Elementformulations
