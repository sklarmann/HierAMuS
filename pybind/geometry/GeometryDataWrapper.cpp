// Copyright 2021 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause

#include "GeometryDataWrapper.h"
#include "geometry/GeometryData.h"
#include "geometry/Volumes.h"
#include "geometry/Special.h"

void HierAMuS::Geometry::GeometryDataWrapper::registerFunctions() {
  this->temp.def(py::init<>())
      .def("requestNewGeometryObject",
           py::overload_cast<HierAMuS::Geometry::GeometryTypes>(
               &HierAMuS::Geometry::GeometryData::requestNewGeometryObject))
      .def("requestNewGeometryObject",
           py::overload_cast<HierAMuS::Geometry::GeometryTypes, indexType>(
               &HierAMuS::Geometry::GeometryData::requestNewGeometryObject))
      .def("getVertex", &HierAMuS::Geometry::GeometryData::getVertex,
           py::return_value_policy::reference)
      .def("getEdge", &HierAMuS::Geometry::GeometryData::getEdge,
           py::return_value_policy::reference)
      .def("getFace", &HierAMuS::Geometry::GeometryData::getFace,
           py::return_value_policy::reference)
      .def("getVolume", &HierAMuS::Geometry::GeometryData::getVolume,
           py::return_value_policy::reference)
      .def("getSpecial", &HierAMuS::Geometry::GeometryData::getSpecial,
           py::return_value_policy::reference)
      .def("print", &HierAMuS::Geometry::GeometryData::print)
      .def("getGeometryElement",
           &HierAMuS::Geometry::GeometryData::getGeometryElement,
           py::return_value_policy::reference)
       .def("checkUpdate",&GeometryData::checkUpdate);
}
