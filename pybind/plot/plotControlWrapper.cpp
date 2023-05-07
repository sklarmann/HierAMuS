// Copyright 2021 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause

#include "plotControlWrapper.h"
#include "plot/plotControl.h"
#include "plot/vtkplotClass.h"
#include "pointercollection/pointercollection.h"
#include "pybind11/stl.h"

void HierAMuS::plotControlWrapper::registerFunctions() {
  this->temp.def(py::init<>())
      .def("initialize", &PlotControl::initialize)
      .def("timeUpdate", &PlotControl::timeUpdate)
      .def("toFile", &PlotControl::toFile);
}
