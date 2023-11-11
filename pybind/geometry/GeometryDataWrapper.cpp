// Copyright 2021 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause

#include "GeometryDataWrapper.h"
#include "geometry/GeometryData.h"
#include "geometry/Volumes/VolumesData.h"
#include "geometry/Special/Special.h"
#include "geometry/VertexData.h"
#include "geometry/Faces/FacesData.h"
#include "geometry/Edges/EdgesData.h"
#include "pointercollection/pointercollection.h"
#include "EquationHandler.h"

void HierAMuS::Geometry::GeometryDataWrapper::registerFunctions() {
  this->temp.def(py::init<>())
      .def("requestNewGeometryObject",
           py::overload_cast<HierAMuS::EquationHandler&,HierAMuS::Geometry::GeometryTypes>(
               &HierAMuS::Geometry::GeometryData::requestNewGeometryObject))
      .def("requestNewGeometryObject",
           py::overload_cast<HierAMuS::EquationHandler&,HierAMuS::Geometry::GeometryTypes, indexType>(
               &HierAMuS::Geometry::GeometryData::requestNewGeometryObject))
      .def("getVertex", &HierAMuS::Geometry::GeometryData::getVertexData,
           py::return_value_policy::reference)
      .def("getEdge", &HierAMuS::Geometry::GeometryData::getEdgeData,
           py::return_value_policy::reference)
      .def("getFace", &HierAMuS::Geometry::GeometryData::getFaceData,
           py::return_value_policy::reference)
      .def("getVolume", &HierAMuS::Geometry::GeometryData::getVolumeData,
           py::return_value_policy::reference)
      .def("getSpecial", &HierAMuS::Geometry::GeometryData::getSpecial,
           py::return_value_policy::reference)
      .def("print", [](GeometryData &self,
                       PointerCollection &pointers) { self.print(pointers.getSPDLogger()); })
      .def("getGeometryElement",
           &HierAMuS::Geometry::GeometryData::getGeometryElement,
           py::return_value_policy::reference)
       .def("checkUpdate",&GeometryData::checkUpdate);
}
