// Copyright 2021 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause

#include "ElementTypesWrapper.h"
#include "finiteElements/NormTypes.h"
#include "pybind/finiteElements/NormTypesWrapper.h"

namespace HierAMuS {
namespace FiniteElement {

void NormtypesWrapper::registerFunctions() {
  this->temp.value("L2Stresses", NormTypes::L2Stresses)
	;
}
} // namespace FiniteElement
}


