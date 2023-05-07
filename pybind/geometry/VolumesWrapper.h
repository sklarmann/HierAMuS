// Copyright 2021 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause

//
// Created by simon on 1/16/22.
//

#pragma once

#include <pybind11/pybind11.h>

namespace py = pybind11;

#include "geometry/Volumes.h"


namespace HierAMuS {
namespace Geometry {
class PyVolumes : public Volumes {
public:
  using Volumes::Volumes;

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
  typedef py::class_<Volumes, PyVolumes, Base, std::shared_ptr<Volumes>>
      pw;
  pw temp;
};
} // namespace Geometry
} // namespace HierAMuS