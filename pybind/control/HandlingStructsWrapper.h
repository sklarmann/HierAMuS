// Copyright 2021 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause

#pragma once
#include <pybind11/pybind11.h>

namespace py = pybind11;

#include "control/HandlingStructs.h"

namespace HierAMuS {
//class PyInfoData : public InfoData {
//public:
//  using InfoData::InfoData;
//
//  // void renew() {
//  //  PYBIND11_OVERRIDE(void, HierAMuS::PointerCollection, renew);
//  //}
//};

class InfoDataWrapper {
public:
  InfoDataWrapper(py::module &m)
      : temp(m, "InfoData"), fileTypes(m, "FileHandling"){};
  void registerFunctions();

private:
  typedef py::class_<InfoData, std::shared_ptr<InfoData>>
      pw;
  pw temp;
  py::enum_<HierAMuS::FileHandling> fileTypes;
};
} // namespace HierAMuS