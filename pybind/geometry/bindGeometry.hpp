// Copyright 2021 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause

#pragma once


#include <pybind11/pybind11.h>

namespace py = pybind11;

void GeometryDataToPybind(py::module &m);
void GeometryTypesToPybind(pybind11::module &m);
void GeomEdgesToPybind(py::module &m);
void GeomLinearEdgeToPybind(py::module &m);
void VertexToPybind(py::module &m);
void BaseToPybind(py::module &m);

void addGeometryFolder(py::module &m) {
  auto GeomModule = m.def_submodule("Geometry", "Geometry Module");
  GeometryTypesToPybind(GeomModule);

  BaseToPybind(GeomModule);
  GeomEdgesToPybind(GeomModule);
  GeomLinearEdgeToPybind(GeomModule);
  VertexToPybind(GeomModule);
  GeometryDataToPybind(GeomModule);
}