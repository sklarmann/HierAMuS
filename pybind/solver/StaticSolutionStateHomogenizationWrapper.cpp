// Copyright 2021 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause

#include "StaticSolutionStateHomogenizationWrapper.h"
#include "pointercollection/pointercollection.h"
#include "control/ParameterList.h"
#include "equations/DegreeOfFreedom.h"
#include "types/MatrixTypes.h"
#include "pybind11/stl.h"
#include "pybind11/eigen.h"


void HierAMuS::StaticSolutionStateHomogenizationWrapper::registerFunctions() {
  this->temp
      .def(py::init<HierAMuS::ParameterList &>())
      .def("setStrains", &StaticSolutionStateHomogenization::setStrains)
      .def("initHomogenization",
           &StaticSolutionStateHomogenization::initHomogenization)
      .def("computeAMatrix", &StaticSolutionStateHomogenization::computeAMatrix)
      .def("homogenize", &StaticSolutionStateHomogenization::homogenize)
      .def("getCMatrix", &StaticSolutionStateHomogenization::getCMatrix)
      .def("getStresses", &StaticSolutionStateHomogenization::getStresses)
      ;
}


