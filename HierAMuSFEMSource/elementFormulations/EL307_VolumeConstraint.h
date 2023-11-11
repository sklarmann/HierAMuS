// Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include "elementFormulations/GenericElementFormulationInterface.h"

namespace HierAMuS {
struct IntegrationPoint;
namespace FiniteElement {
class VolumeConstraint;
}
namespace Geometry {
struct H1Shapes;
}
namespace Elementformulations {

class EL307_VolumeConstraint : public GenericElementFormulationInterface<FiniteElement::VolumeConstraint> {
public:
  EL307_VolumeConstraint(PointerCollection *ptrCol);
  ~EL307_VolumeConstraint() override;
  void readData(PointerCollection &pointers, ParameterList &list) override;
  void setDegreesOfFreedom(PointerCollection &pointers,
                           FiniteElement::VolumeConstraint &elem) override;
  void AdditionalOperations(PointerCollection &pointers,
                            FiniteElement::VolumeConstraint &elem) override;
  auto getDofs(PointerCollection &pointers,
               FiniteElement::GenericFiniteElement *elem)
      -> std::vector<DegreeOfFreedom *> override;
  void setTangentResidual(
      PointerCollection &pointers, FiniteElement::VolumeConstraint &elem,
      Eigen::Matrix<prec, Eigen::Dynamic, Eigen::Dynamic> &stiffness,
      Eigen::Matrix<prec, Eigen::Dynamic, 1> &residual,
      std::vector<DegreeOfFreedom *> &Dofs) override;

  auto getNumberOfIntergrationPoints(PointerCollection &pointers,
                                     FiniteElement::GenericFiniteElement *elem)
      -> indexType override;

  // Paraview
  void toParaviewAdaper(PointerCollection &pointers,
                        FiniteElement::VolumeConstraint &elem,
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


  indexType m_meshIdDisp, m_meshIdLam, m_intOrderDisp, m_mode, m_center;
  prec m_K, m_xmax;

  // prec mu1, mu2;
  // prec e31, e33, e15;
};

} // namespace Elementformulations
} // namespace HierAMuS
