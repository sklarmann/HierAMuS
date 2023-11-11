// Copyright 2021 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause

//
// Created by simon on 1/16/22.
//

#pragma once

#include <pybind11/pybind11.h>

namespace py = pybind11;

#include "geometry/Faces/FacesData.h"


namespace HierAMuS {
namespace Geometry {
class PyFaces : public FacesData {
public:
  using FacesData::FacesData;

  // void renew() {
  //  PYBIND11_OVERRIDE(void, HierAMuS::PointerCollection, renew);
  //}
};

class FacesWrapper {
public:
  FacesWrapper(py::module &m)
      : temp(m, "Faces"){};
  void registerFunctions();

private:
  typedef py::class_<FacesData, PyFaces, GeometryBaseData, std::shared_ptr<FacesData>>
      pw;
  pw temp;
};
} // namespace Geometry
} // namespace HierAMuS