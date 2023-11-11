// Copyright 2021 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause

//
// Created by simon on 1/16/22.
//

#pragma once

#include <pybind11/pybind11.h>

namespace py = pybind11;

#include "geometry/Faces/ScaledBoundary2DData.h"


namespace HierAMuS {
namespace Geometry {
class PyScaledBoundary2D : public ScaledBoundary2DData {
public:
  using ScaledBoundary2DData::ScaledBoundary2DData;

  // void renew() {
  //  PYBIND11_OVERRIDE(void, HierAMuS::PointerCollection, renew);
  //}
};

class ScaledBoundary2DWrapper {
public:
  ScaledBoundary2DWrapper(py::module &m)
      : temp(m, "ScaledBoundary2D"){};
  void registerFunctions();

private:
  typedef py::class_<ScaledBoundary2DData, PyScaledBoundary2D, FacesData,
                     GeometryBaseData, std::shared_ptr<ScaledBoundary2DData>>
      pw;
  pw temp;
};
} // namespace Geometry
} // namespace HierAMuS