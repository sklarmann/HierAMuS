// Copyright 2021 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause

#include "vtkPlotClassWrapper.h"
#include "pointercollection/pointercollection.h"
#include "pybind11/stl.h"

void HierAMuS::vtkPlotClassWrapper::registerFunctions() {
  this->temp.def(py::init<>()).def("timeUpdate",&vtkPlotInterface::timeUpdate)
      .def("initialize",&vtkPlotInterface::initialize)
      .def("toFile", &vtkPlotInterface::toFile)
      .def("addPoint",py::overload_cast<indexType,indexType,indexType,prec,prec,prec>(&vtkPlotInterface::addPoint));
}
