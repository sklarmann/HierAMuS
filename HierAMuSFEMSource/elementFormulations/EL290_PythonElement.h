// Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause


#pragma once


#include "elementFormulations/GenericElementFormulationInterface.h"

#include <tuple>

namespace HierAMuS {
namespace Elementformulations {

class EL290_2DPythonElement : public GenericElementFormulationInterface<FiniteElement::Face> {
public:
  EL290_2DPythonElement(PointerCollection *ptrCol);
  ~EL290_2DPythonElement() override;
  void readData(PointerCollection &pointers, ParameterList &list) override;
  void setDegreesOfFreedom(PointerCollection& pointers, FiniteElement::Face &elem) override;
  void AdditionalOperations(PointerCollection& pointers, FiniteElement::Face &elem) override;
  auto getDofs(PointerCollection& pointers, FiniteElement::GenericFiniteElement *elem)
  -> std::vector<DegreeOfFreedom *> override;
  void setTangentResidual(
    PointerCollection& pointers,
    FiniteElement::Face &elem,
    Eigen::Matrix<prec, Eigen::Dynamic, Eigen::Dynamic> &stiffness,
    Eigen::Matrix<prec, Eigen::Dynamic, 1> &residual, std::vector<DegreeOfFreedom *> &Dofs) override;

  auto getNumberOfIntergrationPoints(PointerCollection &pointers,
                                     FiniteElement::GenericFiniteElement *elem)
      -> indexType override;


  virtual std::tuple<Types::MatrixXX<prec>, Types::VectorX<prec>,
                     std::vector<DegreeOfFreedom *>>
  PySetTangentResidual(FiniteElement::GenericFiniteElement *elem);

  virtual std::tuple<int,int> testfunc(){return {1,1};};

  void setPythonElement(EL290_2DPythonElement *pyElem);

  // plot
  void toParaviewAdaper(PointerCollection &pointers,
                                FiniteElement::Face &elem,
                                vtkPlotInterface &paraviewAdapter,
                                ParaviewSwitch control) override;



private:
  EL290_2DPythonElement *selfPtr;
};

} // namespace Elementformulations
} // namespace HierAMuS
