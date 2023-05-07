// Copyright 2021 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause

#include "StaticSolutionStateWrapper.h"
#include "equations/DegreeOfFreedom.h"
#include "pybind11/stl.h"
#include "pybind11/eigen.h"


void HierAMuS::StaticSolutionStateWrapper::registerFunctions()
{
  this->temp
      .def(py::init<HierAMuS::ParameterList &>())
      .def("setSolver", &HierAMuS::StaticSolutionState::setSolver)
      .def("getSolution",py::overload_cast<indexType>(&StaticSolutionState::getSolution))
      .def("getSolution",py::overload_cast<std::vector<DegreeOfFreedom*>&>(&GenericSolutionState::getSolution))
      .def("energyNorm",&StaticSolutionState::energyNorm)
      .def("residual",&StaticSolutionState::residual)
      ;
}


