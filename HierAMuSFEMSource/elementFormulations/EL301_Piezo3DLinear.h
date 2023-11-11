// Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause


#pragma once

#include "elementFormulations/GenericElementFormulationInterface.h"

namespace HierAMuS::Elementformulations {

class EL301_Piezo3DLinear : public GenericElementFormulationInterface<FiniteElement::Volume> {
public:
  explicit EL301_Piezo3DLinear(PointerCollection *ptrCol);
  ~EL301_Piezo3DLinear() override;
  void readData(PointerCollection &pointers, ParameterList &list) override;
  void setDegreesOfFreedom(PointerCollection &pointers,
                           FiniteElement::Volume &elem) override;
  void AdditionalOperations(PointerCollection &pointers,
                            FiniteElement::Volume &elem) override;
  auto getDofs(PointerCollection& pointers, FiniteElement::GenericFiniteElement *elem) -> std::vector<DegreeOfFreedom*> override;
  void setTangentResidual(
    PointerCollection& pointers,
    FiniteElement::Volume &elem,
    Eigen::Matrix<prec, Eigen::Dynamic, Eigen::Dynamic> &stiffness,
    Eigen::Matrix<prec, Eigen::Dynamic, 1> &residual, std::vector<DegreeOfFreedom *> &Dofs) override;

  auto getNumberOfIntergrationPoints(PointerCollection &pointers,
                                     FiniteElement::GenericFiniteElement *elem)
      -> indexType override;

  // Plot
  void toParaviewAdaper(PointerCollection &pointers,
                        FiniteElement::Volume &elem,
                        vtkPlotInterface &paraviewAdapter,
                        ParaviewSwitch control) override;

private:
  indexType meshIdDisp, meshIdPiezo, intOrderDisp, intOrderPiezo;
};

} // namespace HierAMuS
