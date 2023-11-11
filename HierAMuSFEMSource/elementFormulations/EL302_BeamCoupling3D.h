// Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause


#pragma once

#include "elementFormulations/GenericElementFormulationInterface.h"
#include "shapefunctions/IntegrationsPoints/IntegrationPoints.h"

namespace HierAMuS {
namespace FiniteElement {
class beamInterfaceElement3D;
}
namespace Geometry {
struct H1Shapes;
}
namespace Elementformulations {

class EL302_BeamCoupling3D : public GenericElementFormulationInterface<FiniteElement::beamInterfaceElement3D> {
public:
  EL302_BeamCoupling3D(PointerCollection *ptrCol);
  ~EL302_BeamCoupling3D() override;
  void readData(PointerCollection &pointers, ParameterList &list) override;
  void
  setDegreesOfFreedom(PointerCollection &pointers,
                      FiniteElement::beamInterfaceElement3D &elem) override;
  void AdditionalOperations(PointerCollection& pointers, FiniteElement::beamInterfaceElement3D &elem) override;
  auto getDofs(PointerCollection& pointers, FiniteElement::GenericFiniteElement *elem)
  -> std::vector<DegreeOfFreedom *> override;
  void setTangentResidual(
    PointerCollection& pointers, FiniteElement::beamInterfaceElement3D &elem,
    Eigen::Matrix<prec, Eigen::Dynamic, Eigen::Dynamic> &stiffness,
    Eigen::Matrix<prec, Eigen::Dynamic, 1> &residual, std::vector<DegreeOfFreedom *> &Dofs) override;

  auto getNumberOfIntergrationPoints(PointerCollection& pointers, FiniteElement::GenericFiniteElement *elem)
  -> indexType override;

  // Paraview
  void toParaviewAdaper(PointerCollection &pointers,
                                FiniteElement::beamInterfaceElement3D &elem,
                                vtkPlotInterface &paraviewAdapter,
                                ParaviewSwitch control) override;

private:

  void setTangentResidual_PureDisp(
    PointerCollection& pointers,
    FiniteElement::beamInterfaceElement3D *interfaceElement,
    Eigen::Matrix<prec, Eigen::Dynamic, Eigen::Dynamic> &stiffness,
    Eigen::Matrix<prec, Eigen::Dynamic, 1> &residual, std::vector<DegreeOfFreedom *> &Dofs);


  void setTangentResidual_HuWashizu(
    PointerCollection& pointers,
    FiniteElement::beamInterfaceElement3D *interfaceElement,
    Eigen::Matrix<prec, Eigen::Dynamic, Eigen::Dynamic> &stiffness,
    Eigen::Matrix<prec, Eigen::Dynamic, 1> &residual, std::vector<DegreeOfFreedom *> &Dofs);

  void setTangentResidual_HuWashizu_Warping(
    PointerCollection& pointers,
    FiniteElement::beamInterfaceElement3D *interfaceElement,
    Eigen::Matrix<prec, Eigen::Dynamic, Eigen::Dynamic> &stiffness,
    Eigen::Matrix<prec, Eigen::Dynamic, 1> &residual, std::vector<DegreeOfFreedom *> &Dofs);

    
  void setTangentResidual_Disp_Warping_V1(
    PointerCollection& pointers,
    FiniteElement::beamInterfaceElement3D *interfaceElement,
    Eigen::Matrix<prec, Eigen::Dynamic, Eigen::Dynamic> &stiffness,
    Eigen::Matrix<prec, Eigen::Dynamic, 1> &residual, std::vector<DegreeOfFreedom *> &Dofs);
  


  auto getBMatrixLinear(PointerCollection& pointers, Geometry::H1Shapes &shapes, IntegrationPoint &intPoint, FiniteElement::beamInterfaceElement3D &elem) -> Types::Matrix6X<prec>;


  indexType m_meshIdDisp, m_meshIdRot, m_intOrderDisp, m_mode, m_warpType;
  std::vector<indexType> m_warpBounVerts;

  //prec mu1, mu2;
  //prec e31, e33, e15;
};

} // namespace Elementformulations
} // namespace HierAMuS
