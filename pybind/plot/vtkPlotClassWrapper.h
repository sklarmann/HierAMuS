// Copyright 2021 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause

#pragma once
#include <pybind11/pybind11.h>

namespace py = pybind11;

#include "plot/vtkplotClass.h"

namespace HierAMuS {
//class PyInfoData : public InfoData {
//public:
//  using InfoData::InfoData;
//
//  // void renew() {
//  //  PYBIND11_OVERRIDE(void, HierAMuS::PointerCollection, renew);
//  //}
//};

class vtkPlotClassWrapper {
public:
  vtkPlotClassWrapper(py::module &m) : temp(m, "vtkPlotInterface"){};
  void registerFunctions();

private:
  typedef py::class_<vtkPlotInterface, std::shared_ptr<vtkPlotInterface>>
      pw;
  pw temp;
};
} // namespace HierAMuS