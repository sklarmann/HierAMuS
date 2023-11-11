// Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause



#pragma once


#include <Eigen/Dense>
#include "elementFormulations/GenericElementFormulationInterface.h"

namespace HierAMuS::Elementformulations {

class EL101_Bernoulli2D : public GenericElementFormulationInterface<FiniteElement::Edge> {
public:
  explicit EL101_Bernoulli2D(PointerCollection *ptrCol);
  ~EL101_Bernoulli2D() override;
  void readData(PointerCollection &pointers, ParameterList &list) override;
  void setDegreesOfFreedom(PointerCollection &pointers,
                           FiniteElement::Edge &elem) override;
  void AdditionalOperations(PointerCollection &pointers,
                            FiniteElement::Edge &elem) override;
  auto getDofs(PointerCollection& pointers, FiniteElement::GenericFiniteElement *elem)
  -> std::vector<DegreeOfFreedom *> override;
  void setTangentResidual(
    PointerCollection& pointers,
    FiniteElement::Edge &elem,
    Eigen::Matrix<prec, Eigen::Dynamic, Eigen::Dynamic> &stiffness,
    Eigen::Matrix<prec, Eigen::Dynamic, 1> &residual, std::vector<DegreeOfFreedom *> &Dofs) override;

  auto getNumberOfIntergrationPoints(PointerCollection &pointers,
                                     FiniteElement::GenericFiniteElement *elem)
      -> indexType override;

  // plot
  void toParaviewAdaper(PointerCollection &pointers,
                        FiniteElement::Edge &elem,
                        vtkPlotInterface &paraviewAdapter,
                        ParaviewSwitch control) override;

private:
  indexType meshIdDisp;//, intOrderDisp;
  prec EA, EI;
};
} // namespace HierAMuS
