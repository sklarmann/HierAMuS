// Copyright 2021 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause

#pragma once
#include <pybind11/pybind11.h>

namespace py = pybind11;


#include "TimefunctionWrapper.h"
#include "PropfunctionHandlerWrapper.h"
#include "LoadListWrapper.h"

namespace HierAMuS {
class loadsAddClass {
public:
  loadsAddClass(py::module &m) : func(m), props(m), loads(m)
           {};
  void registerFunctions();

private:
  TimefunctionWrapper func;
  PropfunctionHandlerWrapper props;
  LoadListWrapper loads;
};
} // namespace HierAMuS