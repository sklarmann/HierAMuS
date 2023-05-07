// Copyright 2021 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause

#pragma once
#include <pybind11/pybind11.h>

namespace py = pybind11;

#include "loads/PropfunctionHandler.h"

namespace HierAMuS {
//class PyInfoData : public InfoData {
//public:
//  using InfoData::InfoData;
//
//  // void renew() {
//  //  PYBIND11_OVERRIDE(void, HierAMuS::PointerCollection, renew);
//  //}
//};

class PropfunctionHandlerWrapper {
public:
  PropfunctionHandlerWrapper(py::module &m) : temp(m, "PropfunctionHandler"){};
  void registerFunctions();

private:
  typedef py::class_<PropfunctionHandler, std::shared_ptr<PropfunctionHandler>>
      pw;
  pw temp;
};
} // namespace HierAMuS