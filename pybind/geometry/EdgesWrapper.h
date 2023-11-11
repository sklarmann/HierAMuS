// Copyright 2021 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause

#pragma once
#include <pybind11/pybind11.h>

namespace py = pybind11;

#include "geometry/Edges/EdgesData.h"

namespace HierAMuS {
namespace Geometry {
class PyEdges : public EdgesData {
public:
  using EdgesData::EdgesData;

  // void renew() {
  //  PYBIND11_OVERRIDE(void, HierAMuS::PointerCollection, renew);
  //}
};

class EdgesWrapper {
public:
  EdgesWrapper(py::module &m)
      : temp(m, "Edges"){};
  void registerFunctions();

private:
  typedef py::class_<EdgesData, PyEdges, GeometryBaseData, std::shared_ptr<EdgesData>>
      pw;
  pw temp;
};
} // namespace Geometry
} // namespace HierAMuS