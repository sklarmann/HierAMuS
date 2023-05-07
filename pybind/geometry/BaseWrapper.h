// Copyright 2021 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause

#pragma once
#include <pybind11/pybind11.h>

namespace py = pybind11;

#include "geometry/Base.h"

namespace HierAMuS {
namespace Geometry {
class PyBase : public Base {
public:
  using Base::Base;

  // void renew() {
  //  PYBIND11_OVERRIDE(void, HierAMuS::PointerCollection, renew);
  //}
};

class BaseWrapper {
public:
  BaseWrapper(py::module &m)
      : temp(m, "Base"){};
  void registerFunctions();

private:
  typedef py::class_<Base, PyBase, std::shared_ptr<Base>>
      pw;
  pw temp;
};
} // namespace Geometry
} // namespace HierAMuS