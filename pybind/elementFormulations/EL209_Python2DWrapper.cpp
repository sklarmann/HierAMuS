// Copyright 2021 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause

#include "EL209_Python2DWrapper.h"
#include "pointercollection/pointercollection.h"
#include "control/ParameterList.h"

namespace HierAMuS {
namespace Elementformulations {

void EL290_2DPythonElementWrapper::registerFunctions() {
  this->temp.def(py::init<HierAMuS::PointerCollection *>())
      .def("readData", &EL290_2DPythonElement::readData)
      .def("setDegreesOfFreedom",&EL290_2DPythonElement::setDegreesOfFreedom)
	;
}
} // namespace FiniteElement
} // namespace HierAMuS
