// Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause



#pragma once

#include <forwarddeclaration.h>

#include <Eigen/Dense>
#include <elementFormulations/GenericElementFormulation.h>

namespace HierAMuS::Elementformulations {

class LSFEMBernoulli : public GenericElementFormulation {
public:
  LSFEMBernoulli(PointerCollection *ptrCol);
  ~LSFEMBernoulli() override;
  void readData(PointerCollection &pointers, ParameterList &list) override;
  void setDegreesOfFreedom(PointerCollection& pointers, FiniteElement::GenericFiniteElement *elem) override;
  void setTangentResidual(
    PointerCollection& pointers,
    FiniteElement::GenericFiniteElement *elem,
    Eigen::Matrix<prec, Eigen::Dynamic, Eigen::Dynamic> &stiffness,
    Eigen::Matrix<prec, Eigen::Dynamic, 1> &residual, std::vector<DegreeOfFreedom *> &Dofs) override;
  void setMass(PointerCollection& pointers,
               FiniteElement::GenericFiniteElement *elem,
               Eigen::Matrix<prec, Eigen::Dynamic, Eigen::Dynamic> &stiffness,
               Eigen::Matrix<prec, Eigen::Dynamic, 1> &residual, std::vector<DegreeOfFreedom *> &Dofs) override;
  void getElementsLocalNodalReactions(
      PointerCollection &pointers, FiniteElement::GenericFiniteElement *elem,
      std::map<indexType, std::vector<prec>> &vReacs) override;


  auto getNumberOfIntergrationPoints(PointerCollection &pointers,
                                     FiniteElement::GenericFiniteElement *elem)
      -> indexType override;

  //plot
  void toParaviewAdaper(PointerCollection &pointers,
                        FiniteElement::GenericFiniteElement *elem,
                        vtkPlotInterface &paraviewAdapter,
                        ParaviewSwitch control) override;


private:
  indexType meshIdDisp, propnum, localLoad;
  prec EA, EI, GA, rhoA;
  prec qx, qy, mz;
};

} // namespace HierAMuS
