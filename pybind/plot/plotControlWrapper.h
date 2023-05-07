// Copyright 2021 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause

#pragma once
#include <pybind11/pybind11.h>

namespace py = pybind11;

#include "plot/plotControl.h"

namespace HierAMuS {
//class PyInfoData : public InfoData {
//public:
//  using InfoData::InfoData;
//
//  // void renew() {
//  //  PYBIND11_OVERRIDE(void, HierAMuS::PointerCollection, renew);
//  //}
//};

class plotControlWrapper {
public:
  plotControlWrapper(py::module &m) : temp(m, "PlotControl"){};
  void registerFunctions();

private:
  typedef py::class_<PlotControl, std::shared_ptr<PlotControl>>
      pw;
  pw temp;

};
} // namespace HierAMuS