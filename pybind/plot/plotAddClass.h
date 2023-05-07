// Copyright 2021 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause

#pragma once
#include <pybind11/pybind11.h>

namespace py = pybind11;

#include "vtkPlotClassWrapper.h"
#include "plotControlWrapper.h"

namespace HierAMuS {
class plotAddClass {
public:
  plotAddClass(py::module &m) : vtkplot(m), plotcontrol(m)
           {};
  void registerFunctions();

private:
  vtkPlotClassWrapper vtkplot;
  plotControlWrapper plotcontrol;
};
} // namespace HierAMuS