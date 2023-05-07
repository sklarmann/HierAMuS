// Copyright 2021 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause

#pragma once
#include <pybind11/pybind11.h>

namespace py = pybind11;

#include "loads/LoadList.h"

namespace HierAMuS {
//class PyInfoData : public InfoData {
//public:
//  using InfoData::InfoData;
//
//  // void renew() {
//  //  PYBIND11_OVERRIDE(void, HierAMuS::PointerCollection, renew);
//  //}
//};

class LoadListWrapper {
public:
  LoadListWrapper(py::module &m) : temp(m, "LoadList"){};
  void registerFunctions();

private:
  typedef py::class_<loadList, std::shared_ptr<loadList>>
      pw;
  pw temp;
};
} // namespace HierAMuS