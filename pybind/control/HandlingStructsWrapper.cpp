// Copyright 2021 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause

#include "HandlingStructsWrapper.h"
#include "pybind11/stl.h"

void HierAMuS::InfoDataWrapper::registerFunctions() {
  this->temp.def(py::init<>())
      .def_readonly("Log", &HierAMuS::InfoData::Log)
      .def("setInFile",
           [](InfoData &self, std::string filename) {
             self.fileNames[FileHandling::infile] = filename;
           })
      .def("setOutFile",
           [](InfoData &self, std::string filename) {
             self.fileNames[FileHandling::outfile] = filename;
           })
      .def("setDirectory", [](InfoData &self, std::string dir) {
        self.fileNames[FileHandling::directory] = dir;
      })
      .def("getDir", [](InfoData &self) {
        return self.fileNames[FileHandling::directory];
      });

  this->fileTypes.value("infile", FileHandling::infile)
      .value("outfile", FileHandling::outfile)
      .value("directory", FileHandling::directory);
}
