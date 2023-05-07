// Copyright 2021 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause

#pragma once
#include <pybind11/pybind11.h>

namespace py = pybind11;

#include "NodetypesWrapper.h"
#include "EquationHandlerWrapper.h"
#include "DofStatusWrapper.h"
#include "DegreeOfFreedomWrapper.h"
#include "NodeWrapper.h"

namespace HierAMuS {
class equationsAddClass {
public:
  equationsAddClass(py::module &m) : nodeTypes(m), eqHandler(m), dofStatus(m), Dofs(m), nodes(m)
           {};
  void registerFunctions();

private:
  NodeTypesWrapper nodeTypes;
  EquationHandlerWrapper eqHandler;
  dofStatusWrapper dofStatus;
  DegreeOfFreedomWrapper Dofs;
  GenericNodesWrapper nodes;
};
} // namespace HierAMuS
#undef CNAME