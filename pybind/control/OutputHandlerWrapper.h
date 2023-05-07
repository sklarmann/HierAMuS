// Copyright 2021 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause

#pragma once
#include <pybind11/pybind11.h>

namespace py = pybind11;

#include "control/OutputHandler.h"

namespace HierAMuS {
//class PyOutputHandler : public OutputHandler {
//public:
//  using OutputHandler::OutputHandler;
//
//  // void renew() {
//  //  PYBIND11_OVERRIDE(void, HierAMuS::PointerCollection, renew);
//  //}
//};

class OutputHandlerWrapper {
public:
  OutputHandlerWrapper(py::module &m)
      : temp(m, "OutputHandler"), logl(m, "LogLevel"){};
  void registerFunctions();

private:
  typedef py::class_<OutputHandler, std::shared_ptr<OutputHandler>>
      pw;
  pw temp;
  py::enum_<HierAMuS::LogLevel> logl;
};
} // namespace HierAMuS