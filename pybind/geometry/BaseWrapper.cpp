// Copyright 2021 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause

#include "BaseWrapper.h"
#include "geometry/GeometryBaseData.h"
#include "geometry/GeometryData.h"
#include "pointercollection/pointercollection.h"
#include "LoadList.h"
#include "GenericNodes.h"
#include "pybind11/eigen.h"
#include <pybind11/stl.h>

void HierAMuS::Geometry::BaseWrapper::registerFunctions() {
  this->temp.def("setBoundaryCondition", &GeometryBaseData::setBoundaryCondition)
      .def("setLoad", &GeometryBaseData::setLoad)
      .def("setPrescribedSolution", &GeometryBaseData::setPrescribedSolution)
      .def("getNodes",[](GeometryBaseData &self,PointerCollection &pointers,indexType meshId){
        std::vector<GenericNodes*> nodes;
        self.getNodes(nodes,meshId);
        return nodes;
      },py::return_value_policy::reference)
	;
}
