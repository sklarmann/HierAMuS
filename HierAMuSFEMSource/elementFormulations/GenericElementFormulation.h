// Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause


#pragma once


#include "types/MatrixTypes.h"

#include "finiteElements/NormTypes.h"

#include "solver/HistoryDataNew/HistoryDataStructure.h"

#include <map>
#include <vector>
#include <string>
#include <iostream>


class managementClass;
namespace HierAMuS {
class PointerCollection;
class ParameterList;
class DegreeOfFreedom;
class vtkPlotInterface;
enum class ParaviewSwitch;
class ParameterList;
namespace FiniteElement {
class GenericFiniteElement;
class Edge;
class Face;
class Volume;
class beamInterfaceElement3D;
}

namespace Elementformulations {

class GenericElementFormulation {
public:
  explicit GenericElementFormulation(PointerCollection *ptrCol);
  virtual ~GenericElementFormulation();

  virtual void readData(PointerCollection &pointers, ParameterList &list) = 0;
  virtual void updateRVEHistory(PointerCollection &pointers,
                                FiniteElement::GenericFiniteElement *elem){};
  virtual auto getDofs(PointerCollection &pointers,
                       FiniteElement::GenericFiniteElement *elem)
      -> std::vector<DegreeOfFreedom *> {
    return std::vector<DegreeOfFreedom *>();
  };



  virtual auto computeNorm(PointerCollection &pointers,
                           FiniteElement::GenericFiniteElement *elem,
                           FiniteElement::NormTypes type) -> prec {
    return 1;
  };
  virtual void setMass(PointerCollection &pointers,
                       FiniteElement::GenericFiniteElement *elem,
                       Types::MatrixXX<prec> &stiffness,
                       Types::VectorX<prec> &residual,
                       std::vector<DegreeOfFreedom *> &Dofs){};
  virtual void getElementsLocalNodalReactions(
      PointerCollection &pointers, FiniteElement::GenericFiniteElement *elem,
      std::map<indexType, std::vector<prec>> &vReacs){};

  virtual auto getHistoryDataStructure() -> const HistoryDataStructure &;
  virtual auto
  getNumberOfIntergrationPoints(PointerCollection &pointers,
                                FiniteElement::GenericFiniteElement *elem)
      -> indexType = 0;

  

protected:
  void messageUnprocessed(PointerCollection &pointers, ParameterList &paraMap,
                          std::string elementName);

private:
  const static HistoryDataStructure m_historyDataStructure;
};
} // namespace Elementformulations
} // namespace HierAMuS
