// Copyright 2021 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause

#pragma once
#include <pybind11/pybind11.h>

namespace py = pybind11;

#include "finiteElements/ElementList.h"

namespace HierAMuS {
namespace FiniteElement {
//class PyElementList : public ElementList {
//public:
//  using ElementList::ElementList;
//
//  // void renew() {
//  //  PYBIND11_OVERRIDE(void, HierAMuS::PointerCollection, renew);
//  //}
//};

class ElementListWrapper {
public:
  ElementListWrapper(py::module &m) : temp(m, "ElementList "){};
  void registerFunctions();

private:
  typedef py::class_<ElementList, std::shared_ptr<ElementList>> pw;
  pw temp;
};
} // namespace FiniteElement
} // namespace HierAMuS