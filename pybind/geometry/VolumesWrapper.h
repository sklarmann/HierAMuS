// Copyright 2021 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause

//
// Created by simon on 1/16/22.
//

#pragma once

#include <pybind11/pybind11.h>

namespace py = pybind11;

#include "geometry/Volumes/VolumesData.h"


namespace HierAMuS {
namespace Geometry {
class PyVolumes : public VolumesData {
public:
  using VolumesData::VolumesData;

  // void renew() {
  //  PYBIND11_OVERRIDE(void, HierAMuS::PointerCollection, renew);
  //}
};

class VolumesWrapper {
public:
  VolumesWrapper(py::module &m)
      : temp(m, "Volumes"){};
  void registerFunctions();

private:
  typedef py::class_<VolumesData, PyVolumes, GeometryBaseData,
                     std::shared_ptr<VolumesData>>
      pw;
  pw temp;
};
} // namespace Geometry
} // namespace HierAMuS