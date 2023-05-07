// Copyright 2021 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause

#pragma once
#include <pybind11/pybind11.h>

namespace py = pybind11;

#include "materials/ElementformulationList.h"

namespace HierAMuS {
namespace Materials {

//class PyTransientSolutionNewmark : public MaterialList {
//public:
//  using MaterialList::MaterialList;
//
//  // void renew() {
//  //  PYBIND11_OVERRIDE(void, HierAMuS::PointerCollection, renew);
//  //}
//};

class ElementFormulationListWrapper {
public:
  ElementFormulationListWrapper(py::module &m)
      : temp(m, "ElementFormulationList"){};
  void registerFunctions();

private:
  typedef py::class_<ElementFormulationList, 
                     std::shared_ptr<ElementFormulationList>>
      pw;
  pw temp;
};
} // namespace Materials
} // namespace HierAMuS