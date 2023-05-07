// Copyright 2021 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause

#pragma once
#include <pybind11/pybind11.h>

namespace py = pybind11;

#include "geometry/GeometryData.h"

namespace HierAMuS {
namespace Geometry {
class PyGeometryData : public GeometryData {
public:
  using GeometryData::GeometryData;

  // void renew() {
  //  PYBIND11_OVERRIDE(void, HierAMuS::PointerCollection, renew);
  //}
};

class GeometryDataWrapper {
public:
  GeometryDataWrapper(py::module &m)
      : temp(m, "GeometryData"){};
  void registerFunctions();

private:
  typedef py::class_<GeometryData, PyGeometryData, std::shared_ptr<GeometryData>>
      pw;
  pw temp;
};
} // namespace Geometry
} // namespace HierAMuS