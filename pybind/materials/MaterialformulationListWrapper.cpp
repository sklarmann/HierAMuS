// Copyright 2021 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause

#include "MaterialformulationListWrapper.h"
#include "pointercollection/pointercollection.h"
namespace HierAMuS {
namespace Materials {

void MaterialFormulationListWrapper::registerFunctions() {
  this->temp.def("addMaterial", py::overload_cast<PointerCollection &, indexType,
                                       indexType, ParameterList &>(&MaterialFormulationList::addMaterial))
      .def("addMaterial",
           py::overload_cast<indexType,
                             std::shared_ptr<GenericMaterialFormulation>>(&MaterialFormulationList::addMaterial))
      .def("getMaterial",&MaterialFormulationList::getMaterial)
	;
}
} // namespace Materials
} // namespace HierAMuS
