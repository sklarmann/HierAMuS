// Copyright 2021 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause

#include "GenericSolutionStateWrapper.h"
#include "solver/GenericSolutionState.h"
#include "vector"
#include "PropfunctionHandler.h"
#include "solver/Constraints/ConstraintHandler.h"
#include "control/ParameterList.h"

#include "DegreeOfFreedom.h"

#include "pybind11/stl.h"
#include "pybind11/eigen.h"

void HierAMuS::GenericSolutionStateWrapper::registerFunctions()
{
  this->temp.def(py::init<HierAMuS::ParameterList &>())
      .def("setSparseMatrix",&GenericSolutionState::setSparseMatrix)
      .def("getProps",&GenericSolutionState::getProps)
      .def("nextSolutionStep",&GenericSolutionState::nextSolutionStep)
      .def("assembleSystem", &GenericSolutionState::assembleSystem)
      .def("factorize", &GenericSolutionState::factorize)
      .def("solve", &GenericSolutionState::solve)
      .def("updateSolution", &GenericSolutionState::updateSolution)
      .def("getSolution",py::overload_cast<indexType>(&GenericSolutionState::getSolution))
      .def("getSolution",py::overload_cast<std::vector<DegreeOfFreedom*>&>(&GenericSolutionState::getSolution))
      .def("printSpMat", &GenericSolutionState::printSpMat)
      .def("getNewConstraint", &GenericSolutionState::getNewConstraint)
      .def("getConstraintHandler", &GenericSolutionState::getConstraintHandler,
           py::return_value_policy::reference)
      .def("updateRVEHistory",&GenericSolutionState::updateRVEHistory)
#ifdef USE_SPECTRA
      .def("computeEigenValues",&GenericSolutionState::computeEigenValues)
#endif
      ;
}
