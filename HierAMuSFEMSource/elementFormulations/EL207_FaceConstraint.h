// Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause


#pragma once

#include "shapefunctions/IntegrationsPoints/IntegrationPoints.h"
#include "elementFormulations/GenericElementFormulationInterface.h"


namespace HierAMuS::Geometry{
  struct H1Shapes;
}

namespace HierAMuS::FiniteElement {
class FaceConstraint;
}

namespace HierAMuS::Elementformulations {

class EL207_FaceConstraint : public GenericElementFormulationInterface<FiniteElement::FaceConstraint> {
public:
  explicit EL207_FaceConstraint(PointerCollection *ptrCol);
  ~EL207_FaceConstraint() override;
  void readData(PointerCollection &pointers, ParameterList &list) override;
  void setDegreesOfFreedom(PointerCollection& pointers, FiniteElement::FaceConstraint &elem) override;
  void AdditionalOperations(PointerCollection& pointers, FiniteElement::FaceConstraint &elem) override;
  auto getDofs(PointerCollection& pointers, FiniteElement::GenericFiniteElement *elem)
  -> std::vector<DegreeOfFreedom *> override;
  void setTangentResidual(
    PointerCollection& pointers,
    FiniteElement::FaceConstraint &elem,
    Eigen::Matrix<prec, Eigen::Dynamic, Eigen::Dynamic> &stiffness,
    Eigen::Matrix<prec, Eigen::Dynamic, 1> &residual, std::vector<DegreeOfFreedom *> &Dofs) override;

 
  // plot
  void toParaviewAdaper(PointerCollection &pointers,
                                FiniteElement::FaceConstraint &elem,
                                vtkPlotInterface &paraviewAdapter,
                                ParaviewSwitch control) override;

  auto getHistoryDataStructure() -> const HistoryDataStructure & override;
  auto getNumberOfIntergrationPoints(PointerCollection& pointers, FiniteElement::GenericFiniteElement *elem)
  -> indexType override;
  
private:
  
  void setTangentResidualDispFormulation(
    PointerCollection& pointers,
    FiniteElement::FaceConstraint &elem,
    Eigen::Matrix<prec, Eigen::Dynamic, Eigen::Dynamic> &stiffness,
    Eigen::Matrix<prec, Eigen::Dynamic, 1> &residual, std::vector<DegreeOfFreedom *> &Dofs);
  
  auto getIntegrationPoints(PointerCollection& pointers, FiniteElement::GenericFiniteElement *elem) -> IntegrationPoints;




  indexType meshIdDisp, meshIdRot, meshIdLam, meshIdMu;
  indexType intOrderDisp;
  indexType mode;
  prec k;


  const static HistoryDataStructure m_HistoryDataStructure;
};

} // namespace HierAMuS
