// Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include "MatrixTypes.h"
#include "finiteElements/beamInterfaceElement3D.h"
#include "geometry/Base.h"
#include "pointercollection/pointercollection.h"
#include "shapefunctions/IntegrationsPoints/IntegrationPoints.h"
#include <forwarddeclaration.h>

#include <Eigen/Dense>
#include <elementFormulations/GenericElementFormulation.h>

namespace HierAMuS {
namespace FiniteElement {
class VolumeConstraint;
}
namespace Elementformulations {

class EL307_VolumeConstraint : public GenericElementFormulation {
public:
  EL307_VolumeConstraint(PointerCollection *ptrCol);
  ~EL307_VolumeConstraint() override;
  void readData(PointerCollection &pointers, ParameterList &list) override;
  void setDegreesOfFreedom(PointerCollection &pointers,
                           FiniteElement::GenericFiniteElement *elem) override;
  void AdditionalOperations(PointerCollection &pointers,
                            FiniteElement::GenericFiniteElement *elem) override;
  auto getDofs(PointerCollection &pointers,
               FiniteElement::GenericFiniteElement *elem)
      -> std::vector<DegreeOfFreedom *> override;
  void setTangentResidual(
      PointerCollection &pointers, FiniteElement::GenericFiniteElement *elem,
      Eigen::Matrix<prec, Eigen::Dynamic, Eigen::Dynamic> &stiffness,
      Eigen::Matrix<prec, Eigen::Dynamic, 1> &residual,
      std::vector<DegreeOfFreedom *> &Dofs) override;

  auto getNumberOfIntergrationPoints(PointerCollection &pointers,
                                     FiniteElement::GenericFiniteElement *elem)
      -> indexType override;

  // Paraview
  void toParaviewAdaper(PointerCollection &pointers,
                        FiniteElement::GenericFiniteElement *elem,
                        vtkPlotInterface &paraviewAdapter,
                        ParaviewSwitch control) override;

private:
  void setTangentResidual_XDir(
      PointerCollection &pointers, FiniteElement::VolumeConstraint *volElement,
      Eigen::Matrix<prec, Eigen::Dynamic, Eigen::Dynamic> &stiffness,
      Eigen::Matrix<prec, Eigen::Dynamic, 1> &residual,
      std::vector<DegreeOfFreedom *> &Dofs);

  void setTangentResidual_TempGradient(
      PointerCollection &pointers, FiniteElement::VolumeConstraint *volElement,
      Eigen::Matrix<prec, Eigen::Dynamic, Eigen::Dynamic> &stiffness,
      Eigen::Matrix<prec, Eigen::Dynamic, 1> &residual,
      std::vector<DegreeOfFreedom *> &Dofs);

  auto getBMatrixLinear(PointerCollection &pointers, Geometry::H1Shapes &shapes,
                        IntegrationPoint &intPoint,
                        FiniteElement::beamInterfaceElement3D &elem)
      -> Types::Matrix6X<prec>;

  indexType m_meshIdDisp, m_meshIdLam, m_intOrderDisp, m_mode, m_center;
  prec m_K, m_xmax;

  // prec mu1, mu2;
  // prec e31, e33, e15;
};

} // namespace Elementformulations
} // namespace HierAMuS
