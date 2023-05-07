// Copyright 2021 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause

#pragma once
#include <pybind11/pybind11.h>

namespace py = pybind11;

#include "MaterialListWrapper.h"
#include "MaterialWrapper.h"
#include "ElementformulationListWrapper.h"
#include "MaterialformulationListWrapper.h"
#include "GenericMaterialFormulationWrapper.h"

namespace HierAMuS {
class materialAddClass {
public:
  materialAddClass(py::module &m)
      : matlist(m), material(m), elemList(m), matFormList(m), GenericMat(m)
	{};
  void registerFunctions();

private:
  Materials::MaterialListWrapper matlist;
  Materials::MaterialWrapper material;
  Materials::ElementFormulationListWrapper elemList;
  Materials::MaterialFormulationListWrapper matFormList;
  Materials::GenericMaterialFormulationWrapper GenericMat;
};
} // namespace HierAMuS
#undef CNAME