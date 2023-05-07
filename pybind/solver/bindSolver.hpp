// Copyright 2021 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <pybind11/pybind11.h>

namespace py = pybind11;

void GenericSolutionStateToPybind(pybind11::module &m);
void SolutionTypesToPybind(pybind11::module &m);
void StaticSolutionStateToPybind(pybind11::module &m);
void SolverTypesToPybind(pybind11::module &m);

void addSolverFolder(py::module &m) {
  SolutionTypesToPybind(m);
  SolverTypesToPybind(m);
  GenericSolutionStateToPybind(m);
  StaticSolutionStateToPybind(m);
}