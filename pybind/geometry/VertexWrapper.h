// Copyright 2021 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause

#pragma once
#include <pybind11/pybind11.h>

namespace py = pybind11;

#include "geometry/Vertex.h"
#include "geometry/Base.h"

namespace HierAMuS {
namespace Geometry {
class PyVertex : public Vertex {
public:
  using Vertex::Vertex;

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
  typedef py::class_<Vertex, PyVertex, Base, std::shared_ptr<Vertex>>
      pw;
  pw temp;
};
} // namespace Geometry
} // namespace HierAMuS