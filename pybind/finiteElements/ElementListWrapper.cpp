// Copyright 2021 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause

#include "ElementListWrapper.h"

namespace HierAMuS {
namespace FiniteElement {

void ElementListWrapper::registerFunctions() {
  this->temp.def(py::init<>())
      .def("requestNewElement", &ElementList::requestNewElement,py::return_value_policy::reference)
      .def("getElement", &ElementList::getElement,py::return_value_policy::reference)
      .def("getNumberOfElements", &ElementList::getNumberOfElements)
      .def("setDegreesOfFreedom", &ElementList::setDegreesOfFreedom);
}
} // namespace FiniteElement
} // namespace HierAMuS
