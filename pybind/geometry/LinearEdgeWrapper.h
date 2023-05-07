// Copyright 2021 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause

#pragma once
#include <pybind11/pybind11.h>

namespace py = pybind11;

#include "geometry/LinearEdge.h"

namespace HierAMuS {
namespace Geometry {
class PyLinearEdge : public LinearEdge {
public:
  using LinearEdge::LinearEdge;

  // void renew() {
  //  PYBIND11_OVERRIDE(void, HierAMuS::PointerCollection, renew);
  //}
};

class LinearEdgeWrapper {
public:
  LinearEdgeWrapper(py::module &m)
      : temp(m, "LinearEdge "){};
  void registerFunctions();

private:
  typedef py::class_<LinearEdge, PyLinearEdge, Edges, Base, std::shared_ptr<LinearEdge>>
      pw;
  pw temp;
};
} // namespace Geometry
} // namespace HierAMuS