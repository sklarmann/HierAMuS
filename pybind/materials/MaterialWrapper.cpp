// Copyright 2021 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause

#include "MaterialWrapper.h"
#include "elementFormulations/GenericElementFormulation.h"
#include "materials/GenericMaterialFormulation.h"
#include "pointercollection/pointercollection.h"

namespace HierAMuS {
namespace Materials {

void MaterialWrapper::registerFunctions() {
  this->temp.def("setElementForumaltion", &Material::setElementForumaltion)
      .def("setMaterialFormulation",&Material::setMaterialFormulation)
      .def("getElementFormulation",&Material::getElementFormulation)
      .def("getMaterialFormulation",&Material::getMaterialFormulation)
	;
}
} // namespace Materials
} // namespace HierAMuS
