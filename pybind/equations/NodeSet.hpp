// Copyright 2021 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause

#pragma once
#include <memory>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;

#include "equations/NodeSet.h"

// class PyEquationHandler : public HierAMuS::EquationHandler {
// public:
//  using HierAMuS::EquationHandler::EquationHandler;
//
//  //void renew() {
//  //  PYBIND11_OVERRIDE(void, HierAMuS::PointerCollection, renew);
//  //}
//};

void NodeSetToPybind(pybind11::module &m) {
  py::class_<HierAMuS::NodeSet,
             std::shared_ptr<HierAMuS::NodeSet>>(
      m, "NodeSet")
      .def(py::init<>())
      .def("setType",&HierAMuS::NodeSet::setType)
      .def("setMeshId",&HierAMuS::NodeSet::setMeshId)
      .def("getMeshId",&HierAMuS::NodeSet::getMeshId)
      .def("setNumberOfNodes",&HierAMuS::NodeSet::setNumberOfNodes)
      .def("getNumberOfNodes",&HierAMuS::NodeSet::getNumberOfNodes)
      .def("getNodeSetType",&HierAMuS::NodeSet::getNodeSetType)
      .def("getNodes",&HierAMuS::NodeSet::getNodes)
      .def("getDegreesOfFreedom",&HierAMuS::NodeSet::getDegreesOfFreedom)
	;
}
