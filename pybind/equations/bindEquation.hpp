// Copyright 2021 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause

#pragma once
#include "EquationHandler.hpp"
#include "DofStatus.hpp"
#include "NodeSet.hpp"
#include "DegreeOfFreedom.hpp"

#include <pybind11/pybind11.h>


namespace py = pybind11;

void addEquationFolder(py::module &m) {
  EquationHandlerToPybind(m);
  DofStatusToPybind(m);
  NodeSetToPybind(m);
  DegreeOfFreedomToPybind(m);
}