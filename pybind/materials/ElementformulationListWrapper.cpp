// Copyright 2021 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause

#include "ElementformulationListWrapper.h"
#include "pointercollection/pointercollection.h"
#include "control/ParameterList.h"

namespace HierAMuS {
namespace Materials {

void ElementFormulationListWrapper::registerFunctions() {
  this->temp.def("addElementFormulation",
                 py::overload_cast<PointerCollection &, indexType, indexType,
                                   ParameterList&>(
                     &ElementFormulationList::addElementFormulation))
      .def("addElementFormulation", py::overload_cast<
                   indexType,
                   std::shared_ptr<
                       Elementformulations::GenericElementFormulation>>(&ElementFormulationList::addElementFormulation))
      .def("getElementFormulation",&ElementFormulationList::getElementFormulation)
	;
}
} // namespace Materials
} // namespace HierAMuS
