// Copyright 2021 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause

#pragma once
#include <pybind11/pybind11.h>

namespace py = pybind11;

#include "HandlingStructsWrapper.h"
#include "OutputHandlerWrapper.h"
#include "ParameterListWrapper.h"



namespace HierAMuS {
class controlAddClass {
public:
  controlAddClass(py::module &m) : infoData(m), outputData(m), paramList(m)
           {};
  void registerFunctions();

private:
  InfoDataWrapper infoData;
  OutputHandlerWrapper outputData;
  ParameterListWrapper paramList;
};
} // namespace HierAMuS