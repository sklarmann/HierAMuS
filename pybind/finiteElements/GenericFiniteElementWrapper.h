// Copyright 2021 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause

#pragma once
#include <pybind11/pybind11.h>

namespace py = pybind11;

#include "finiteElements/GenericFiniteElement.h"

namespace HierAMuS {
namespace FiniteElement {
class PyGenericFiniteElement : public GenericFiniteElement {
public:
  using GenericFiniteElement::GenericFiniteElement;

  // void renew() {
  //  PYBIND11_OVERRIDE(void, HierAMuS::PointerCollection, renew);
  //}
};

class GenericFiniteElementWrapper {
public:
  GenericFiniteElementWrapper(py::module &m) : temp(m, "GenericFiniteElement"){};
  void registerFunctions();

private:
  typedef py::class_<GenericFiniteElement,
                     std::shared_ptr<GenericFiniteElement>>
      pw;
  pw temp;
};
} // namespace FiniteElement
} // namespace HierAMuS