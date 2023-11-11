// Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include "elementFormulations/GenericElementFormulationInterface.h"
#include "finiteElements/GenericFiniteElement.h"
#include "pointercollection/pointercollection.h"
#include "materials/Material.h"

#include "control/OutputHandler.h"

namespace HierAMuS::FiniteElement {

template <class T>
class GenericFiniteElementInterface : public GenericFiniteElement {

public:
  GenericFiniteElementInterface() = default;
  virtual ~GenericFiniteElementInterface() = default;

  void GenericSetDegreesOfFreedom(PointerCollection &pointers) override {

    this->getElementFormulation(pointers)->setDegreesOfFreedom(
        pointers, static_cast<T &>(*this));
  }

  void GenericAdditionalOperations(PointerCollection &pointers) override {
    this->getElementFormulation(pointers)->AdditionalOperations(
        pointers, static_cast<T &>(*this));
  }

  auto getDofs(PointerCollection &pointers)
      -> std::vector<DegreeOfFreedom *> override {
    return this->m_material->getElementFormulation(pointers)->getDofs(
        pointers, static_cast<T *>(this));
  }

  void
  GenericSetTangentResidual(PointerCollection &pointers,
                            Types::MatrixXX<prec> &stiffness,
                            Types::VectorX<prec> &residual,
                            std::vector<DegreeOfFreedom *> &Dofs) override {

    this->getElementFormulation(pointers)->setTangentResidual(
        pointers, static_cast<T &>(*this), stiffness, residual, Dofs);
  }

  void toParaviewAdapter(PointerCollection &pointers,
                         vtkPlotInterface &catalyst,
                         ParaviewSwitch ToDo) override {
    this->getElementFormulation(pointers)->toParaviewAdaper(
        pointers, static_cast<T &>(*this), catalyst, ToDo);
  }

private:
  auto getElementFormulation(PointerCollection &pointers)
      -> Elementformulations::GenericElementFormulationInterface<T> * {
    auto eformPtr = this->m_material->getElementFormulation(pointers);

    auto retPtr = dynamic_cast<
        Elementformulations::GenericElementFormulationInterface<T> *>(
        &(*eformPtr));

    if (retPtr != nullptr)
      return retPtr;
    pointers.getSPDLogger().critical(
        "Error when trying to cast Elementformulation. Associated element "
        "formulation cannot handle the specific finite element!");
    throw std::runtime_error(
        "Mismatch between element and elementformulation!");
  }
};
} // namespace HierAMuS::FiniteElement
