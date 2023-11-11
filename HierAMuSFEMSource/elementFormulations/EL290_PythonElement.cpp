// Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause

#include "elementFormulations/EL290_PythonElement.h"

#include "finiteElements/Face.h"

namespace HierAMuS {
namespace Elementformulations {
EL290_2DPythonElement::EL290_2DPythonElement(PointerCollection *ptr)
    : GenericElementFormulationInterface(ptr) {}

EL290_2DPythonElement::~EL290_2DPythonElement() {}

void EL290_2DPythonElement::readData(PointerCollection &pointers,
                                     ParameterList &list) {}

void EL290_2DPythonElement::setDegreesOfFreedom(PointerCollection &pointers,
                                                FiniteElement::Face &elem) {}

void EL290_2DPythonElement::AdditionalOperations(
    PointerCollection& pointers, FiniteElement::Face &elem) {}

auto EL290_2DPythonElement::getDofs(PointerCollection &pointers,
                                    FiniteElement::GenericFiniteElement *elem)
    -> std::vector<DegreeOfFreedom *> {
  return {};
}

void EL290_2DPythonElement::setTangentResidual(
    PointerCollection &pointers, FiniteElement::Face &elem,
    Types::MatrixXX<prec> &stiffness, Types::VectorX<prec> &residual,
    std::vector<DegreeOfFreedom *> &Dofs) {

  auto [s, r, d] = this->PySetTangentResidual(&elem);
  stiffness = s;
  residual = r;
  Dofs = d;
}

auto EL290_2DPythonElement::getNumberOfIntergrationPoints(
    PointerCollection &pointers, FiniteElement::GenericFiniteElement *elem)
    -> indexType {
  return 0;
}

std::tuple<Types::MatrixXX<prec>, Types::VectorX<prec>,
           std::vector<DegreeOfFreedom *>>
EL290_2DPythonElement::PySetTangentResidual(
    FiniteElement::GenericFiniteElement *elem) {

  return std::tuple<Types::MatrixXX<prec>, Types::VectorX<prec>,
                    std::vector<DegreeOfFreedom *>>();
}

void EL290_2DPythonElement::setPythonElement(EL290_2DPythonElement *pyElem) {
  this->selfPtr = pyElem;
}

void EL290_2DPythonElement::toParaviewAdaper(
    PointerCollection &pointers, FiniteElement::Face &elem,
    vtkPlotInterface &paraviewAdapter, ParaviewSwitch control) {}

} // namespace Elementformulations
} // namespace HierAMuS