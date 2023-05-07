// Copyright 2021 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause

#include "BaseWrapper.h"
#include "geometry/Base.h"
#include "geometry/GeometryData.h"
#include "pointercollection/pointercollection.h"
#include "equations/GenericNodes.h"
#include "pybind11/eigen.h"
#include <pybind11/stl.h>

void HierAMuS::Geometry::BaseWrapper::registerFunctions() {
  this->temp.def("setBoundaryCondition", &Base::setBoundaryCondition)
      .def("setLoad", &Base::setLoad)
      .def("setPrescribedSolution", &Base::setPrescribedSolution)
      .def("getNodes",[](Base &self,PointerCollection &pointers,indexType meshId){
        std::vector<GenericNodes*> nodes;
        self.getNodes(pointers,nodes,meshId);
        return nodes;
      },py::return_value_policy::reference)
      .def("getCoordinates",py::overload_cast<>(&Base::getCoordinates))
      .def("getVertNums",[](Base &self){
        std::vector<indexType> vertsOut;
        self.getVerts(vertsOut);
        return vertsOut;
      })
	;
}
