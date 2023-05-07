// Copyright 2021 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause

#pragma once
#include <pybind11/pybind11.h>

namespace py = pybind11;

#include "ElementListWrapper.h"
#include "ElementTypesWrapper.h"
#include "GenericFiniteElementWrapper.h"
#include "NormTypesWrapper.h"

namespace HierAMuS {
class finiteElementsAddClass {
public:
  finiteElementsAddClass(py::module &m)
      : elemList(m), genericElem(m), elemTypes(m), normTypes(m){};
  void registerFunctions();

private:
  FiniteElement::ElementListWrapper elemList;
  FiniteElement::GenericFiniteElementWrapper genericElem;
  FiniteElement::ElementtypesWrapper elemTypes;
  FiniteElement::NormtypesWrapper normTypes;
};
} // namespace HierAMuS
#undef CNAME