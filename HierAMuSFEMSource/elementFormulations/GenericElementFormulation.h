// Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause


#pragma once

#include <forwarddeclaration.h>

#include "plot/vtkplotClassBase.h"

//#include <Eigen/Dense>
#include <control/ParameterList.h>
#include <map>
#include <vector>

#include "finiteElements/NormTypes.h"

#include "solver/HistoryDataNew/HistoryDataStructure.h"

class managementClass;
namespace HierAMuS::Elementformulations {

class GenericElementFormulation {
public:
  explicit GenericElementFormulation(PointerCollection *ptrCol);
  virtual ~GenericElementFormulation();

  virtual void readData(PointerCollection &pointers, ParameterList &list);
  virtual void setDegreesOfFreedom(PointerCollection& pointers, FiniteElement::GenericFiniteElement *elem);
  virtual void AdditionalOperations(PointerCollection& pointers, FiniteElement::GenericFiniteElement *elem);
  virtual void updateRVEHistory(PointerCollection &pointers,
                                FiniteElement::GenericFiniteElement *elem){};
  virtual auto getDofs(PointerCollection& pointers, FiniteElement::GenericFiniteElement *elem)
  -> std::vector<DegreeOfFreedom *> {
    return std::vector<DegreeOfFreedom *>();
  };
  virtual void setTangentResidual(
    PointerCollection& pointers,
    FiniteElement::GenericFiniteElement *elem,
    Eigen::Matrix<prec, Eigen::Dynamic, Eigen::Dynamic> &stiffness,
    Eigen::Matrix<prec, Eigen::Dynamic, 1> &residual, std::vector<DegreeOfFreedom *> &Dofs){};

  virtual auto computeNorm(PointerCollection &pointers,
                           FiniteElement::GenericFiniteElement *elem,
                           FiniteElement::NormTypes type) -> prec {
    return 1;
  };
  virtual void
  setMass(PointerCollection& pointers,
          FiniteElement::GenericFiniteElement *elem,
          Eigen::Matrix<prec, Eigen::Dynamic, Eigen::Dynamic> &stiffness,
          Eigen::Matrix<prec, Eigen::Dynamic, 1> &residual, std::vector<DegreeOfFreedom *> &Dofs){};
  virtual void getElementsLocalNodalReactions(
      PointerCollection &pointers, FiniteElement::GenericFiniteElement *elem,
      std::map<indexType, std::vector<prec>> &vReacs){};

  virtual auto getHistoryDataStructure() -> const HistoryDataStructure &;
  virtual auto
  getNumberOfIntergrationPoints(PointerCollection& pointers, FiniteElement::GenericFiniteElement *elem)
  -> indexType = 0;

  // plot
  virtual void toParaviewAdaper(PointerCollection &pointers,
                                FiniteElement::GenericFiniteElement *elem,
                                vtkPlotInterface &paraviewAdapter,
                                ParaviewSwitch control) = 0;


protected:
  void messageUnprocessed(PointerCollection &pointers, ParameterList &paraMap, std::string elementName);

private:
  const static HistoryDataStructure m_historyDataStructure;
};
} // namespace HierAMuS::Elementformulations
