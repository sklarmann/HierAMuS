// Copyright 2021 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause

#pragma once
#include <pybind11/pybind11.h>

namespace py = pybind11;

#include "geometry/Edges/LinearEdgeData.h"

namespace HierAMuS {
namespace Geometry {
class PyLinearEdge : public LinearEdgeData {
public:
  using LinearEdgeData::LinearEdgeData;

  // void renew() {
  //  PYBIND11_OVERRIDE(void, HierAMuS::PointerCollection, renew);
  //}
};

class LinearEdgeWrapper {
public:
  LinearEdgeWrapper(py::module &m)
      : temp(m, "LinearEdge"){};
  void registerFunctions();

private:
  typedef py::class_<LinearEdgeData, PyLinearEdge, EdgesData, GeometryBaseData,
                     std::shared_ptr<LinearEdgeData>>
      pw;
  pw temp;
};
} // namespace Geometry
} // namespace HierAMuS