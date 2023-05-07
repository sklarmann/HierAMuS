// Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause



#pragma once

#include <forwarddeclaration.h>

#include <elementFormulations/GenericElementFormulation.h>
//#include <Eigen/Dense>

namespace HierAMuS::Elementformulations {

class EL202_Piezo2D : public GenericElementFormulation {
public:
  explicit EL202_Piezo2D(PointerCollection *ptrCol);
  ~EL202_Piezo2D() override;
  void readData(PointerCollection &pointers, ParameterList &list) override;
  void setDegreesOfFreedom(PointerCollection& pointers, FiniteElement::GenericFiniteElement *elem) override;
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

  //Plot
  void toParaviewAdaper(PointerCollection &pointers,
                        FiniteElement::GenericFiniteElement *elem,
                        vtkPlotInterface &paraviewAdapter,
                        ParaviewSwitch control) override;

private:
  indexType meshIdDisp, intOrderDisp;
  prec emod, nu;
  prec mu1, mu2;
  prec e31, e33, e15;
};

} // namespace HierAMuS
