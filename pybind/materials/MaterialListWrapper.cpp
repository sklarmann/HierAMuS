// Copyright 2021 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause

#include "MaterialListWrapper.h"
#include "materials/Material.h"

namespace HierAMuS {
namespace Materials {

void MaterialListWrapper::registerFunctions() {
  this->temp.def("getMaterial", &MaterialList::getMaterial,py::return_value_policy::reference)
	.def("getNumberOfMaterials",&MaterialList::getNumberOfMaterials);
}
} // namespace Materials
} // namespace HierAMuS
