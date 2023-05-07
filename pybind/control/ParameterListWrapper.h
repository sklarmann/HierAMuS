// Copyright 2021 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause

#pragma once
#include <pybind11/pybind11.h>

namespace py = pybind11;

#include "control/ParameterList.h"

namespace HierAMuS {
//class PyParameterList : public ParameterList {
//public:
//  using ParameterList::ParameterList;
//
//  // void renew() {
//  //  PYBIND11_OVERRIDE(void, HierAMuS::PointerCollection, renew);
//  //}
//};

class ParameterListWrapper {
public:
  ParameterListWrapper(py::module &m)
      : temp(m, "ParameterList"){};
  void registerFunctions();

private:
  typedef py::class_<ParameterList, std::shared_ptr<ParameterList>>
      pw;
  pw temp;
};
} // namespace HierAMuS