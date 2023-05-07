// Copyright 2021 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause

#include "OutputHandlerWrapper.h"
#include "control/OutputHandler.h"

void HierAMuS::OutputHandlerWrapper::registerFunctions() {
  this->temp.def(py::init<>())
      .def("setLogLevel", &HierAMuS::OutputHandler::setLogLevel)
      .def("openLogFile", &HierAMuS::OutputHandler::openLogFile)
      .def("closeLogFile", &HierAMuS::OutputHandler::closeLogFile)
      .def("outputString", &HierAMuS::OutputHandler::outputString)
      .def("outputStringBasic", &HierAMuS::OutputHandler::outputStringBasic)
      .def("outputStringDebug", &HierAMuS::OutputHandler::outputStringDebug)
      ;

	this->logl.value("NoLog", HierAMuS::LogLevel::NoLog)
      .value("BasicLog", HierAMuS::LogLevel::BasicLog)
      .value("FullLog", HierAMuS::LogLevel::FullLog);

}
