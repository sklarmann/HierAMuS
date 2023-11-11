// Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause


#pragma once

#include "elementFormulations/GenericElementFormulationInterface.h"
#include "shapefunctions/IntegrationsPoints/IntegrationPoints.h"

namespace HierAMuS {
namespace Geometry {
struct H1Shapes;
}
namespace Elementformulations {

class EL304_QPVolumeElement : public GenericElementFormulationInterface<FiniteElement::Volume> {
public:
  EL304_QPVolumeElement(PointerCollection *ptrCol);
  ~EL304_QPVolumeElement() override;
  void readData(PointerCollection &pointers, ParameterList &list) override;
  void setDegreesOfFreedom(PointerCollection& pointers, FiniteElement::Volume &elem) override;
  void AdditionalOperations(PointerCollection& pointers, FiniteElement::Volume &elem) override;
  auto getDofs(PointerCollection& pointers, FiniteElement::GenericFiniteElement *elem) -> std::vector<DegreeOfFreedom*> override;
  void setTangentResidual(
    PointerCollection& pointers,
    FiniteElement::Volume &elem,
    Eigen::Matrix<prec, Eigen::Dynamic, Eigen::Dynamic> &stiffness,
    Eigen::Matrix<prec, Eigen::Dynamic, 1> &residual, std::vector<DegreeOfFreedom *> &Dofs) override;

  // Paraview
  void toParaviewAdaper(PointerCollection &pointers,
                                FiniteElement::Volume &elem,
                                vtkPlotInterface &paraviewAdapter,
                                ParaviewSwitch control) override;
  auto getHistoryDataStructure() -> const HistoryDataStructure & override;
  auto getNumberOfIntergrationPoints(PointerCollection& pointers, FiniteElement::GenericFiniteElement *elem)
  -> indexType override;
  

private:
  void setTangentResidualNonLinear(
    PointerCollection& pointers,
    FiniteElement::Volume &elem,
    Eigen::Matrix<prec, Eigen::Dynamic, Eigen::Dynamic> &stiffness,
    Eigen::Matrix<prec, Eigen::Dynamic, 1> &residual, std::vector<DegreeOfFreedom *> &Dofs);

  auto getIntegrationPoints(PointerCollection &pointers,
    FiniteElement::GenericFiniteElement& elem) -> IntegrationPoints;

  void getBdBv(PointerCollection &pointers, Types::Matrix6X<prec> &Bd,
               Types::VectorXT<prec> &Bv,
               Geometry::H1Shapes &shapes);

  const static HistoryDataStructure m_HistoryDataStructure;

  indexType meshIdDisp, dispOrder, pressureOrder; 
  prec mu, kappa;
};

} // namespace Elementformulations
} // namespace HierAMuS
