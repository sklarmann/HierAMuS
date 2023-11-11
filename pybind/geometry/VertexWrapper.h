// Copyright 2021 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause

#pragma once
#include <pybind11/pybind11.h>

namespace py = pybind11;

#include "geometry/VertexData.h"
#include "geometry/GeometryBaseData.h"

namespace HierAMuS {
namespace Geometry {
class PyVertex : public VertexData {
public:
  using VertexData::VertexData;

  // void renew() {
  //  PYBIND11_OVERRIDE(void, HierAMuS::PointerCollection, renew);
  //}
};

class VertexWrapper {
public:
  VertexWrapper(py::module &m)
      : temp(m, "Vertex "){};
  void registerFunctions();

private:
  typedef py::class_<VertexData, PyVertex, GeometryBaseData,
                     std::shared_ptr<VertexData>>
      pw;
  pw temp;
};
} // namespace Geometry
} // namespace HierAMuS