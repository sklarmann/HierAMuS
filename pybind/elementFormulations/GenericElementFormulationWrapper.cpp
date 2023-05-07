// Copyright 2021 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause

#include "GenericElementFormulationWrapper.h"
#include "pointercollection/pointercollection.h"

namespace HierAMuS {
namespace Elementformulations {

void GenericElementFormulationWrapper::registerFunctions() {
  this->temp//.def(py::init<HierAMuS::PointerCollection *>())
      .def("readData", &GenericElementFormulation::readData)
	;
}
} // namespace FiniteElement
} // namespace HierAMuS
