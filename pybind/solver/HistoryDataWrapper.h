// Copyright 2021 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause

#pragma once
#include <pybind11/pybind11.h>

namespace py = pybind11;

#include "solver/HistoryData.h"

namespace HierAMuS {
//class PyHistoryData : public HierAMuS::HistoryData {
//public:
//  using HierAMuS::HistoryData::HistoryData;
//
//  // void renew() {
//  //  PYBIND11_OVERRIDE(void, HierAMuS::PointerCollection, renew);
//  //}
//};

class HistoryDataWrapper {
public:
  HistoryDataWrapper(py::module &m)
      : temp(m, "HistoryData"){};
  void registerFunctions();

private:
  typedef py::class_<HistoryData, 
                     std::shared_ptr<HistoryData>>
      pw;
  pw temp;
};
} // namespace HierAMuS