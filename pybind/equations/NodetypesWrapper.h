// Copyright 2021 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause

#pragma once
#include <pybind11/pybind11.h>

namespace py = pybind11;

#include "equations/Nodetypes.h"

namespace HierAMuS {

class NodeTypesWrapper {
public:
  NodeTypesWrapper(py::module &m) : temp(m, "NodeTypes"){};
  void registerFunctions();

private:
  typedef py::enum_<NodeTypes> pw;
  pw temp;
};
} // namespace HierAMuS