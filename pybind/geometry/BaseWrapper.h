// Copyright 2021 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause

#pragma once
#include <pybind11/pybind11.h>

namespace py = pybind11;

#include "geometry/GeometryBaseData.h"

namespace HierAMuS {
namespace Geometry {
class PyGeometryBaseData : public GeometryBaseData {
public:
  using GeometryBaseData::GeometryBaseData;

  // void renew() {
  //  PYBIND11_OVERRIDE(void, HierAMuS::PointerCollection, renew);
  //}
};

class BaseWrapper {
public:
  BaseWrapper(py::module &m)
      : temp(m, "GeometryBaseData"){};
  void registerFunctions();

private:
  typedef py::class_<GeometryBaseData, PyGeometryBaseData,
                     std::shared_ptr<GeometryBaseData>>
      pw;
  pw temp;
};
} // namespace Geometry
} // namespace HierAMuS