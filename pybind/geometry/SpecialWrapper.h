// Copyright 2021 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause

//
// Created by simon on 1/16/22.
//

#pragma once

#include <pybind11/pybind11.h>

namespace py = pybind11;

#include "geometry/Special/Special.h"


namespace HierAMuS {
namespace Geometry {
class PySpecial : public Special {
public:
  using Special::Special;

  // void renew() {
  //  PYBIND11_OVERRIDE(void, HierAMuS::PointerCollection, renew);
  //}
};

class SpecialWrapper {
public:
  SpecialWrapper(py::module &m)
      : temp(m, "Special"){};
  void registerFunctions();

private:
  typedef py::class_<Special, PySpecial, GeometryBaseData,
                     std::shared_ptr<Special>>
      pw;
  pw temp;
};
} // namespace Geometry
} // namespace HierAMuS