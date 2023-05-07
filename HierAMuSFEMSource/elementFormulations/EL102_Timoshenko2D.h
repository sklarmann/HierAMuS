// Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause



#pragma once

#include "types/MatrixTypes.h"
#include <forwarddeclaration.h>

#include <elementFormulations/GenericElementFormulation.h>
#include <Eigen/Dense>

namespace HierAMuS::Elementformulations {

class EL102_Timoshenko2D : public GenericElementFormulation {
public:
  explicit EL102_Timoshenko2D(PointerCollection *ptrCol);
  ~EL102_Timoshenko2D() override;
  void readData(PointerCollection &pointers, ParameterList &list) override;
  void setDegreesOfFreedom(PointerCollection& pointers, FiniteElement::GenericFiniteElement *elem) override;
  void AdditionalOperations(PointerCollection& pointers, FiniteElement::GenericFiniteElement *elem) override;
  auto getDofs(PointerCollection& pointers, FiniteElement::GenericFiniteElement *elem)
  -> std::vector<DegreeOfFreedom *> override;
  void setTangentResidual(
    PointerCollection& pointers,
    FiniteElement::GenericFiniteElement *elem,
    Eigen::Matrix<prec, Eigen::Dynamic, Eigen::Dynamic> &stiffness,
    Eigen::Matrix<prec, Eigen::Dynamic, 1> &residual, std::vector<DegreeOfFreedom *> &Dofs) override;
  void setTangentResidualLinear(
    PointerCollection& pointers,
    FiniteElement::GenericFiniteElement *elem,
    Eigen::Matrix<prec, Eigen::Dynamic, Eigen::Dynamic> &stiffness,
    Eigen::Matrix<prec, Eigen::Dynamic, 1> &residual, std::vector<DegreeOfFreedom *> &Dofs);
  void setTangentResidualNonlinear(
    PointerCollection& pointers,
    FiniteElement::GenericFiniteElement *elem,
    Eigen::Matrix<prec, Eigen::Dynamic, Eigen::Dynamic> &stiffness,
    Eigen::Matrix<prec, Eigen::Dynamic, 1> &residual, std::vector<DegreeOfFreedom *> &Dofs);
  void setMass(PointerCollection& pointers,
               FiniteElement::GenericFiniteElement *elem,
               Eigen::Matrix<prec, Eigen::Dynamic, Eigen::Dynamic> &stiffness,
               Eigen::Matrix<prec, Eigen::Dynamic, 1> &residual, std::vector<DegreeOfFreedom *> &Dofs) override;
  void getElementsLocalNodalReactions(
      PointerCollection &ptrCol, FiniteElement::GenericFiniteElement *elem,
      std::map<indexType, std::vector<prec>> &vReacs) override;

  // plot
  void toParaviewAdaper(PointerCollection &pointers,
                                FiniteElement::GenericFiniteElement *elem,
                                vtkPlotInterface &paraviewAdapter,
                                ParaviewSwitch control) override;
  

  auto getHistoryDataStructure() -> const HistoryDataStructure & override;
  auto getNumberOfIntergrationPoints(PointerCollection& pointers, FiniteElement::GenericFiniteElement *elem)
  -> indexType override;
  

private:
  Types::Matrix33<prec> getMaterialMatrix() const;
  indexType meshIdDisp, meshIdRot, disporder, rotorder, propnum, localLoad;
  indexType mode;
  prec EA, EI, GA, rhoA;
  prec qx, qy, mz;

  
  const static HistoryDataStructure m_HistoryDataStructure;
};

} // namespace HierAMuS
