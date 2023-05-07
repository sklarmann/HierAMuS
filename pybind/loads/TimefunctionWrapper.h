// Copyright 2021 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause

#pragma once
#include <pybind11/pybind11.h>

namespace py = pybind11;

#include "loads/Timefunction.h"

namespace HierAMuS {
//class PyInfoData : public InfoData {
//public:
//  using InfoData::InfoData;
//
//  // void renew() {
//  //  PYBIND11_OVERRIDE(void, HierAMuS::PointerCollection, renew);
//  //}
//};

class TimefunctionWrapper {
public:
  TimefunctionWrapper(py::module &m) : temp(m, "Timefunction"), fun(m, "Function"){};
  void registerFunctions();

private:
  typedef py::class_<Timefunction, std::shared_ptr<Timefunction>>
      pw;
  pw temp;
  py::class_<Function, std::shared_ptr<Function>> fun;
};
} // namespace HierAMuS